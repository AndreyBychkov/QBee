import math
import configparser
import pickle
import numpy as np
from sympy.polys.rings import PolyElement
from sympy import *
from typing import Callable, List, Optional, Set, Collection
from functools import partial
from .selection import *  # replace with .selection if you want pip install
from .util import *  # replace with .util if you want pip install

from memory_profiler import profile

config = configparser.ConfigParser({
    'logging_enable': False,
    'progress_bar_enable': False,
    'logging_file': 'log/log.feather',
    'quad_systems_file': 'log/quad_systems.pkl'
})
config.read("../config.ini")
log_enable = eval(config.get('DEFAULT', 'logging_enable'))  # Security Error here, but does not matter I believe
pb_enable = eval(config.get('DEFAULT', 'progress_bar_enable'))  # Security Error here, but does not matter I believe
log_file = config.get('DEFAULT', 'logging_file')
quad_systems_file = config.get('DEFAULT', 'quad_systems_file')
if log_enable:
    print(f"Log file will be produced as {log_file}, quadratizations will be saved as {quad_systems_file}")


def quadratize(polynomials: List[PolyElement],
               selection_strategy=aeqd_strategy,
               pruning_functions: Optional[Union[Tuple, List]] = None):
    if pruning_functions is None:
        pruning_functions = (pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound)
    system = PolynomialSystem(polynomials)
    algo = BranchAndBound(system, selection_strategy, (pruning_by_best_nvars,) + pruning_functions)
    quad_res = algo.quadratize()
    if pb_enable:
        print("Quadratized system:", quad_res)
    quad_system = apply_quadratization(polynomials, quad_res.system.introduced_vars)
    return quad_system


# ------------------------------------------------------------------------------

class PolynomialSystem:
    def __init__(self, polynomials: List[PolyElement]):
        """
        polynomials - right-hand sides of the ODE system listed in the same order as
                      the variables in the polynomial ring
        """
        self.dim = len(polynomials[0].ring.gens)
        self.gen_syms = list(map(lambda g: sp.Symbol(str(g)), polynomials[0].ring.gens))

        # put not monomials but differences in the exponents between rhs and lhs
        self.rhs = dict()
        for i, p in enumerate(polynomials):
            self.rhs[i] = set()
            for m in p.to_dict().keys():
                mlist = list(m)
                mlist[i] -= 1
                self.rhs[i].add(tuple(mlist))

        self.vars, self.squares, self.nonsquares = set(), set(), set()
        self.add_var(tuple([0] * self.dim))
        for i in range(self.dim):
            self.add_var(tuple([1 if i == j else 0 for j in range(self.dim)]))
        self.original_degree = max(map(monomial_deg, self.nonsquares))

    @property
    def introduced_vars(self):
        return tuple(filter(lambda v: monomial_deg(v) >= 2, self.vars))

    def add_var(self, v):
        for i in range(self.dim):
            if v[i] > 0:
                for m in self.rhs[i]:
                    self.nonsquares.add(tuple([v[j] + m[j] for j in range(self.dim)]))

        self.vars.add(v)
        for u in self.vars:
            self.squares.add(tuple([u[i] + v[i] for i in range(self.dim)]))
        self.nonsquares = set(filter(lambda s: s not in self.squares, self.nonsquares))

    def copy(self):
        return pickle.loads(pickle.dumps(self, -1))

    def is_quadratized(self):
        return not self.nonsquares

    def get_smallest_nonsquare(self):
        return min([(np.prod([d + 1 for d in m]), m) for m in self.nonsquares])[1]

    def next_generation(self, strategy=default_strategy):
        if len(self.nonsquares) == 0:
            return list()
        new_gen = []
        for d in get_decompositions(self.get_smallest_nonsquare()):
            c = pickle.loads(pickle.dumps(self, -1))  # inline self.copy for speedup
            for v in d:
                c.add_var(v)
            new_gen.append(c)

        return sorted(new_gen, key=strategy)

    def new_vars_count(self):
        return len(self.vars) - self.dim - 1

    def __str__(self):
        return f"{self._introduced_variables_str()}"

    def __repr__(self):
        return f"{self._introduced_variables_str()}"

    def _introduced_variables_str(self):
        """Representation for visualization"""
        return sorted(map(partial(monom2str, gens=self.gen_syms), self.introduced_vars))


# ------------------------------------------------------------------------------

class QuadratizationResult:
    def __init__(self,
                 system: Optional[PolynomialSystem],
                 introduced_vars: int,
                 nodes_traversed: int):
        self.system = system
        self.introduced_vars = introduced_vars
        self.nodes_traversed = nodes_traversed

    def __repr__(self):
        if self.system is None:
            return "No quadratization found under the given condition\n" + \
                   f"Nodes traversed: {self.nodes_traversed}"
        return f"Number of introduced variables: {self.introduced_vars}\n" + \
               f"Nodes traversed: {self.nodes_traversed}\n" + \
               f"Introduced variables: {self.system}"


# ------------------------------------------------------------------------------

Pruning = Callable[..., bool]


class Algorithm:
    def __init__(self,
                 poly_system: PolynomialSystem,
                 strategy: SelectionStrategy = default_strategy,
                 pruning_funcs: Collection[Pruning] = None):
        self._system = poly_system
        self._strategy = strategy
        self._pruning_funs = list(pruning_funcs) if pruning_funcs is not None else [lambda a, b, *_: False]
        self._nodes_traversed = 0
        self.preliminary_upper_bound = math.inf

        if len(list(poly_system.vars)[0]) > 1:
            points = list()
            for i, eq in poly_system.rhs.items():
                for m in eq:
                    new_m = list(m)
                    new_m[i] = new_m[i] if new_m[i] > 0 else 0
                    points.append(tuple(new_m))
            self.dominating_monomials = set()
            for m in points + list(poly_system.vars):
                if not dominated(m, self.dominating_monomials):
                    self.dominating_monomials.add(m)

    def quadratize(self) -> QuadratizationResult:
        pass

    def traverse_all(self, to_depth: int, pred: Callable[[PolynomialSystem], bool]):
        res = set()
        self._dls(self._system, to_depth, pred, res)
        self._final_iter()
        return res

    @progress_bar(is_stop=False, enabled=pb_enable)
    def _dls(self, part_res: PolynomialSystem, to_depth: int, pred: Callable[[PolynomialSystem], bool], res: set):
        if part_res.new_vars_count() > to_depth:
            return

        if pred(part_res):
            res.add(part_res)
        else:
            for next_system in self.next_gen(part_res):
                self._dls(next_system, to_depth, pred, res)
        return

    @dump_results(log_enable, quad_systems_file)
    def get_optimal_quadratizations(self) -> Set[PolynomialSystem]:
        optimal_first = self.quadratize()
        print(optimal_first)
        return self.get_quadratizations(optimal_first.introduced_vars)

    @dump_results(log_enable, quad_systems_file)
    def get_quadratizations(self, depth: int) -> Set[PolynomialSystem]:
        return self.traverse_all(depth, lambda s: s.is_quadratized())

    def add_pruning(self, termination_criteria: Pruning) -> None:
        self._pruning_funs.append(termination_criteria)

    @property
    def selection_strategy(self):
        return self._strategy

    @selection_strategy.setter
    def selection_strategy(self, value):
        self._strategy = value

    @logged(log_enable, log_file)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.selection_strategy)

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        pass


# ------------------------------------------------------------------------------

class BranchAndBound(Algorithm):

    def domination_upper_bound(self):
        system = self._system.copy()
        algo = BranchAndBound(system, aeqd_strategy,
                              [pruning_by_best_nvars,
                               pruning_by_quadratic_upper_bound,
                               pruning_by_squarefree_graphs,
                               partial(pruning_by_domination, dominators=self.dominating_monomials)])
        res = algo.quadratize()
        upper_bound = res.introduced_vars
        self.preliminary_upper_bound = upper_bound
        if upper_bound != math.inf:
            self.add_pruning(partial(pruning_by_vars_number, nvars=upper_bound))

    @timed(enabled=pb_enable)
    def quadratize(self, cond: Callable[[PolynomialSystem], bool] = lambda _: True) -> QuadratizationResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, self.preliminary_upper_bound, cond)
        self._final_iter()
        self._save_results(opt_system)
        return QuadratizationResult(opt_system, nvars, traversed)

    @progress_bar(is_stop=False, enabled=pb_enable)
    def _bnb_step(self, part_res: PolynomialSystem, best_nvars, cond) \
            -> Tuple[Union[int, float], Optional[PolynomialSystem], int]:
        self._nodes_traversed += 1
        if part_res.is_quadratized() and cond(part_res):
            return part_res.new_vars_count(), part_res, 1
        if any(map(lambda f: f(self, part_res, best_nvars), self._pruning_funs)):
            return math.inf, None, 1

        traversed_total = 1
        min_nvars, best_system = best_nvars, None
        for next_system in self.next_gen(part_res):
            nvars, opt_system, traversed = self._bnb_step(next_system, min_nvars, cond)
            traversed_total += traversed
            if nvars < min_nvars:
                min_nvars = nvars
                best_system = opt_system
        return min_nvars, best_system, traversed_total

    @logged(log_enable, log_file)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.selection_strategy)

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @dump_results(enabled=True, log_file=quad_systems_file)
    def _save_results(self, opt_system):
        return [opt_system, ]


# ------------------------------------------------------------------------------

def pruning_by_nodes_processed(algo: Algorithm, _: PolynomialSystem, *args, nodes_processed: int):
    """
    Stops a search if it's computer 'nodes_processed' nodes.

    :examples
        >>> from functools import partial
        >>> pruning = partial(pruning_by_nodes_processed, nodes_processed=100000)
    """
    if algo._nodes_traversed >= nodes_processed:
        return True
    return False


def pruning_by_elapsed_time(algo: Algorithm, system: PolynomialSystem, *args, start_t, max_t):
    """
    Stops a search if 'max_t' was exceeded.

    :examples
        >>> from functools import partial
        >>> pruning = partial(pruning_by_elapsed_time, start_t=time(), max_t=100) # 100 seconds
    :return:
    """
    curr_t = time()
    if curr_t - start_t >= max_t:
        return True
    return False


def pruning_by_vars_number(_: Algorithm, system: PolynomialSystem, *args, nvars: int):
    """
    Search for quadratization with at most 'nvars' components

    :Examples
        >>> from functools import partial
        >>> pruning = partial(pruning_by_vars_number, nvars=10)

    """
    if system.new_vars_count() >= nvars:
        return True
    return False


def pruning_by_best_nvars(a: Algorithm, part_res: PolynomialSystem, *args):
    """Branch-and-Bound default pruning """
    best_nvars, *_ = args
    if part_res.new_vars_count() >= best_nvars - 1:
        return True
    return False


def pruning_by_quadratic_upper_bound(a: Algorithm, part_res: PolynomialSystem, *args):
    best_nvars, *_ = args

    degree_one_monomials = dict()
    for ns in part_res.nonsquares:
        for v in part_res.vars:
            diff = tuple([ns[i] - v[i] for i in range(part_res.dim)])
            if not any([x < 0 for x in diff]):
                if diff in degree_one_monomials:
                    degree_one_monomials[diff] += 1
                else:
                    degree_one_monomials[diff] = 1
    degree_one_count = sorted(degree_one_monomials.values(), reverse=True)

    needed_new_vars = 0
    degree_two_monoms = len(part_res.nonsquares)
    while True:
        if needed_new_vars < len(degree_one_count):
            degree_two_monoms -= degree_one_count[needed_new_vars]
        needed_new_vars += 1
        if degree_two_monoms <= (needed_new_vars * (needed_new_vars + 1)) // 2:
            break

    lower_bound = part_res.new_vars_count() + needed_new_vars

    if lower_bound >= best_nvars:
        return True
    return False


# the (i, j)-th element is the maximal number of edges in a graph
# on i vertices with at most j loops without a four-cycle
# the first items in each row can be checked against Thm 2 from https://doi.org/10.1002/jgt.3190130107
MAX_C4_FREE_EDGES = [
    [0],
    [0, 1],
    [1, 2, 2],
    [3, 3, 4, 4],
    [4, 5, 5, 6, 6],
    [6, 6, 7, 7, 8, 8],
    [7, 8, 9, 9, 9, 10, 10],
    [9, 10, 11, 12, 12, 12, 12, 12]
]


def pruning_by_squarefree_graphs(a: Algorithm, part_res: PolynomialSystem, *args):
    best_nvars, *_ = args

    no_C4_monoms = set()
    sums_of_monoms = set()
    ns_ordered = sorted([m for m in part_res.nonsquares if sum(m) % 2 == 1], key=sum, reverse=True) + \
                 sorted([m for m in part_res.nonsquares if sum(m) % 2 == 0], key=sum, reverse=True)
    for m in ns_ordered:
        new_sums = set()
        to_add = True
        for mm in no_C4_monoms.union(set([m])):
            s = tuple([m[i] + mm[i] for i in range(len(m))])
            if s in sums_of_monoms:
                to_add = False
                break
            new_sums.add(s)
        if to_add:
            sums_of_monoms = sums_of_monoms.union(new_sums)
            no_C4_monoms.add(m)

    no_C4_edges = len(no_C4_monoms)
    max_loops = sum([1 for m in no_C4_monoms if all([i % 2 == 0 for i in m])])

    degree_one_monomials = dict()
    for ns in no_C4_monoms:
        for v in part_res.vars:
            diff = tuple([ns[i] - v[i] for i in range(part_res.dim)])
            if not any([x < 0 for x in diff]):
                if diff in degree_one_monomials:
                    degree_one_monomials[diff] += 1
                else:
                    degree_one_monomials[diff] = 1
    degree_one_count = sorted(degree_one_monomials.values(), reverse=True)

    needed_new_vars = 0
    while True:
        if needed_new_vars < len(degree_one_count):
            no_C4_edges -= degree_one_count[needed_new_vars]
        needed_new_vars += 1
        if needed_new_vars >= len(MAX_C4_FREE_EDGES):
            break
        if no_C4_edges <= MAX_C4_FREE_EDGES[needed_new_vars][min(max_loops, needed_new_vars)]:
            break

    lower_bound = part_res.new_vars_count() + needed_new_vars

    if lower_bound >= best_nvars:
        return True

    return False


def pruning_by_domination(a: BranchAndBound, part_res: PolynomialSystem, *args, dominators):
    if not part_res.introduced_vars:
        return False

    return not all([dominated(m, dominators) for m in part_res.introduced_vars])


# ------------------------------------------------------------------------------------------------

def with_higher_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) > system.original_degree, system.vars))


def with_le_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) <= system.original_degree, system.vars))
