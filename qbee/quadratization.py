from __future__ import annotations

import copy
import math
import configparser
import signal
import pickle
import numpy as np
from sympy.polys.rings import PolyElement
from sympy.core.function import AppliedUndef
from queue import Queue
from typing import Callable, List, Optional, Set, Collection
from ordered_set import OrderedSet
from functools import partial
from operator import add
from .selection import *  # replace with .selection if you want pip install
from .util import *  # replace with .util if you want pip install
from .polynomialization import EquationSystem, polynomialize
from .printer import print_qbee, str_qbee

from memory_profiler import profile

config = configparser.ConfigParser({
    'logging_enable': False,
    'progress_bar_enable': True,
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
               conditions: Collection["SystemCondition"] = (),
               calc_upper_bound=True,
               generation_strategy=default_generation,
               scoring: Scoring = default_scoring,
               pruning_functions: Collection["Pruning"] | None = None, new_vars_name='w',
               start_new_vars_with=0) -> QuadratizationResult | None:
    """
    Quadratize a system of ODEs with the polynomial right-hand side.

    :param conditions:
    :param polynomials: List of polynomials that are the right-hand side of a system and are built from the elements of sympy.PolyRing.
     Left-hand side is given according to the definition of variables in sympy.ring.
    :param selection_strategy: UPDATE
    :param pruning_functions: predicates that remove transformations from the search space
    :param new_vars_name: base name for new variables. Example: new_var_name='z' => z0, z1, z2, ...
    :param start_new_vars_with: Initial index for new variables. Example: start_new_vars_with=3 => w3, w4, ...
    :return: quadratized system or None if there is none found

    Example:
        >>> from sympy import ring, QQ
        >>> R, x, y = ring("x, y", QQ)
        >>> quad_res = quadratize([x**2 * y, x * y**3],new_vars_name='z',start_new_vars_with=1)
        >>> print(quad_res)
        ==================================================
        Quadratization result
        ==================================================
        Number of introduced variables: 2
        Nodes traversed: 16
        Introduced variables:
        z{1} = x*y**2
        z{2} = x*y
        x' = x*z{2}
        y' = y*z{1}
        z{1}' = 2*z{1}**2 + z{1}*z{2}
        z{2}' = z{1}*z{2} + z{2}**2

    """
    if pruning_functions is None:
        pruning_functions = default_pruning_rules
    system = PolynomialSystem(polynomials)
    algo = BranchAndBound(system, conditions, generation_strategy, scoring, (pruning_by_best_nvars,) + tuple(pruning_functions))
    if calc_upper_bound:
        algo.domination_upper_bound()
    algo_res = algo.quadratize()
    if pb_enable:
        print("=" * 50)
        print("Quadratization result")
        print("=" * 50)
        print(algo_res.print(new_vars_name, start_new_vars_with))
        print()
    if algo_res.system is not None:
        quad_eqs, eq_vars = apply_quadratization(polynomials, algo_res.system.introduced_vars,
                                                 new_vars_name, start_new_vars_with)
        return QuadratizationResult(quad_eqs, eq_vars, algo_res)
    return None


def polynomialize_and_quadratize_ode(system: Union[EquationSystem, List[Tuple[sp.Symbol, sp.Expr]]],
                                     input_der_orders=None, conditions: Collection["SystemCondition"] = (),
                                     polynomialization_upper_bound=10, calc_upper_bound=True,
                                     generation_strategy=default_generation,
                                     scoring: Scoring = default_scoring,
                                     pruning_functions: Collection["Pruning"] | None = None,
                                     new_vars_name="w_", start_new_vars_with=0) -> Optional[QuadratizationResult]:

    """
    Polynomialize and then quadratize a system of ODEs with the continuous right-hand side.

    :param system: system of equations in the form [(X, f(X)), ...] where the left-hand side is the derivatives.
    :param input_der_orders: mapping of input variables to maximum order of their derivatives. For example {T: 2} => T in C2
    :param new_vars_name: base name for new variables. Example: new_var_name='z' => z0, z1, z2, ...
    :param start_new_vars_with: initial index for new variables. Example: start_new_vars_with=3 => w3, w4, ...
    :return: quadratized system or None if there is none found

    Example:
        >>> from qbee import *
        >>> from sympy import exp
        >>> x, y, u = functions("x, y, u")
        >>> p = parameters("p")
        >>> quad_res = polynomialize_and_quadratize_ode([(x, y / (1 + exp(-p * x))), (y, x * exp(y) + u)],input_der_orders={u: 0},new_vars_name='z',start_new_vars_with=1)
        >>> print(quad_res)
        Variables introduced in polynomialization:
        z{1} = exp(-p*x)
        z{2} = 1/(z{1} + 1)
        z{3} = exp(y)
        ==================================================
        Quadratization result
        ==================================================
        Number of introduced variables: 4
        Nodes traversed: 116
        Introduced variables:
        z{4} = y*z{2}
        z{5} = y*z{1}*z{2}**2
        z{6} = z{1}*z{2}**2
        z{7} = x*z{3}
        x' = z{4}
        y' = u + z{7}
        z{1}' = -p*z{1}*z{4}
        z{2}' = p*z{4}*z{6}
        z{3}' = u*z{3} + z{3}*z{7}
        u' = 0
        z{4}' = p*z{4}*z{5} + u*z{2} + z{2}*z{7}
        z{5}' = -p*z{4}*z{5} + 2*p*z{5}**2 + u*z{6} + z{6}*z{7}
        z{6}' = -p*z{4}*z{6} + 2*p*z{5}*z{6}
        z{7}' = u*z{7} + z{3}*z{4} + z{7}**2
    """
    if input_der_orders is None:
        input_der_orders = dict()
    if pb_enable:
        # TODO: temporary solution, should incorporate printing variables with non-integer powers into subs. equations
        print("Variables introduced in polynomialization:")
    poly_system = polynomialize(system, polynomialization_upper_bound,
                                new_var_name=new_vars_name, start_new_vars_with=start_new_vars_with)
    if pb_enable:
        poly_system.print_substitutions()
    poly_equations, excl_inputs = poly_system.to_poly_equations(input_der_orders)
    without_excl_inputs = partial(without_variables, excl_vars=excl_inputs)
    pruning_by_decl_inputs = partial(pruning_by_declining_variables, excl_vars=excl_inputs)
    if pruning_functions is None:
        pruning_functions = default_pruning_rules
    quad_result = quadratize(poly_equations,
                             conditions=[without_excl_inputs, *conditions],
                             calc_upper_bound=calc_upper_bound,
                             generation_strategy=partial(generation_strategy, excl_vars=excl_inputs),
                             scoring=scoring,
                             pruning_functions=[pruning_by_best_nvars, pruning_by_decl_inputs, *pruning_functions],
                             new_vars_name=new_vars_name,
                             start_new_vars_with=start_new_vars_with + len(poly_system) - len(system))
    if quad_result:
        quad_result.polynomialization = poly_system
    return quad_result


def polynomialize_and_quadratize(start_system: List[Tuple[sp.Symbol, sp.Expr]], input_der_orders: Optional[Dict] = None,
                                 conditions: Collection["SystemCondition"] = (), polynomialization_upper_bound=10,
                                calc_quadr_upper_bound=True,
                                 generation_strategy=default_generation,
                                 scoring: Scoring = default_scoring,
                                 pruning_functions: Collection["Pruning"] | None = None,
                                 new_vars_name="w_", start_new_vars_with=0) -> Optional[QuadratizationResult]:
    queue = Queue()
    queue.put(start_system)
    if input_der_orders is None:
        inputs = select_inputs(start_system)
        input_der_orders = {i: 0 for i in inputs}
    while not queue.empty():
        system = queue.get_nowait()
        inputs_pde = select_pde_inputs(system)
        input_orders_with_pde = {i: 0 for i in inputs_pde}
        input_orders_with_pde.update(input_der_orders)
        if pb_enable:
            print("Current spatial time derivatives equations:")
            print("...")
            for eq in system[len(start_system):]:
                print(f"{str_qbee(eq[0])} = {str_qbee(eq[1])}")
            print()

        quad_res = polynomialize_and_quadratize_ode(system, input_orders_with_pde, conditions,
                                                    polynomialization_upper_bound,
                                                    calc_upper_bound=calc_quadr_upper_bound,
                                                    generation_strategy=generation_strategy, scoring=scoring,
                                                    pruning_functions=pruning_functions, new_vars_name=new_vars_name,
                                                    start_new_vars_with=start_new_vars_with)
        if quad_res:
            return quad_res
        for i in inputs_pde:
            new_sys = deepcopy(system)
            ex, dx = rm_last_diff(i)
            try:
                new_sys.append((i, get_rhs(system, ex).diff(dx)))
                queue.put(new_sys)
            except AttributeError as e:
                pass
    return None


def get_rhs(system, sym):
    for lhs, rhs in system:
        if lhs == sym:
            return rhs
    return None


def rm_last_diff(der):
    if len(der.variables) == 1:
        return der.expr, der.variables[0]
    elif len(der.variables) > 1:
        return sp.Derivative(der.expr, *der.variables[:-1]), der.variables[-1]


def select_inputs(system):
    lhs, rhs = zip(*system)
    funcs = set(reduce(lambda l, r: l | r, [eq.atoms(AppliedUndef) for eq in rhs]))
    lhs_args = set(sp.flatten([eq.args for eq in lhs if not isinstance(eq, sp.Derivative)]))
    return set(filter(lambda f: (f not in lhs) and (f not in lhs_args), funcs))


def select_pde_inputs(system):
    lhs, rhs = zip(*system)
    return set(filter(lambda v: v not in lhs, reduce(lambda l, r: l | r, [eq.atoms(sp.Derivative) for eq in rhs])))


# ------------------------------------------------------------------------------

class PolynomialSystem:
    def __init__(self, polynomials: List[PolyElement]):
        """
        polynomials - right-hand sides of the ODE system listed in the same order as
                      the variables in the polynomial ring
        """
        gens = polynomials[0].ring.gens
        self.dim = len(gens)
        self.gen_symbols = list(map(lambda g: sp.Symbol(str(g)), gens))

        # put not monomials but differences in the exponents between rhs and lhs
        self.rhs = dict()
        self.laurent = False
        for i, p in enumerate(polynomials):
            self.rhs[i] = set()
            for m in p.to_dict().keys():
                mlist = list(m)
                if any([x < 0 for x in mlist]):
                    self.laurent = True
                mlist[i] -= 1
                self.rhs[i].add(tuple(mlist))

        # PERFORMANCE DEGRADATION: Ordered set is used for correct indexing in output,
        # but the performance degradation is 5-10%
        self.vars, self.squares, self.nonsquares = OrderedSet(), set(), set()
        self.add_var(tuple([0] * self.dim))
        for i in range(self.dim):
            self.add_var(tuple([1 if i == j else 0 for j in range(self.dim)]))
        self.original_degree = max(map(monomial_deg, self.nonsquares.union(self.squares)))

    @property
    def introduced_vars(self):
        return tuple(filter(lambda v: sum(map(abs, v)) >= 2 or sum(v) < 0, self.vars))

    def add_var(self, v):
        for i in range(self.dim):
            if v[i] != 0:
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
        return min([(np.prod([abs(d) + 1 for d in m]), m) for m in self.nonsquares])[1]

    def next_generation(self, generation=default_generation, scoring=default_scoring):
        if len(self.nonsquares) == 0:
            return list()
        new_gen = []
        for d in generation(self):
            c = pickle.loads(pickle.dumps(self, -1))  # inline self.copy for speedup
            for v in d:
                c.add_var(v)
            new_gen.append(c)

        return sorted(new_gen, key=scoring)

    def new_vars_count(self):
        return len(self.vars) - self.dim - 1

    def to_str(self, new_var_name='z_', start_id=0):
        return '\n'.join([
            new_var_name + ("{%d}" % i) + " = " + monom2str(m, self.gen_symbols)
            for i, m in enumerate(self.introduced_vars, start_id)
        ])

    def __str__(self):
        return f"{self._introduced_variables_str()}"

    def __repr__(self):
        return f"{self._introduced_variables_str()}"

    def _introduced_variables_str(self):
        """Representation for visualization"""
        return sorted(map(partial(monom2str, gens=self.gen_symbols), self.introduced_vars))


# ------------------------------------------------------------------------------

class AlgorithmResult:
    def __init__(self,
                 system: Optional[PolynomialSystem],
                 introduced_vars: int,
                 nodes_traversed: int):
        self.system = system
        self.num_introduced_vars = introduced_vars
        self.nodes_traversed = nodes_traversed

    def print(self, new_var_name="z_", start_new_vars_with=0):
        if self.system is None:
            return "No quadratization found under the given condition\n" + \
                   f"Nodes traversed: {self.nodes_traversed}"
        return f"Number of introduced variables: {self.num_introduced_vars}\n" + \
               f"Nodes traversed: {self.nodes_traversed}\n" + \
               "Introduced variables:\n" + self.system.to_str(new_var_name, start_new_vars_with)

    def __repr__(self):
        return self.print()


class QuadratizationResult:
    def __init__(self, equations, variables, quad_res: AlgorithmResult, poly_res: EquationSystem | None = None):
        self.nodes_traversed = quad_res.nodes_traversed
        self.rhs = [copy.copy(e) for e in equations]
        self.lhs = derivatives(variables)
        self.quadratization = quad_res.system
        self.polynomialization = poly_res

    def to_list(self):
        return [self[i] for i in range(len(self.rhs))]

    def introduced_variables_str(self):
        if self.polynomialization:
            base_name = self.polynomialization.variables.base_var_name
            quad_start_index = self.polynomialization.variables.start_new_vars_with + \
                               len(self.polynomialization.variables.generated)
            return self.polynomialization.substitution_equations_str() + \
                   '\n' + \
                   self.quadratization.to_str(base_name, quad_start_index)
        else:
            return self.quadratization.to_str()

    def __getitem__(self, i):
        return sp.Eq(self.lhs[i], self.rhs[i])

    def __repr__(self):
        return '\n'.join([
            f"{dx} = {fx}" for dx, fx in zip(self.lhs, self.rhs)
        ])

    @property
    def new_vars_count(self):
        return self.quadratization.new_vars_count() +\
               (len(self.polynomialization.variables.generated) if self.polynomialization else 0)

    def __len__(self):
        return len(self.lhs)


# ------------------------------------------------------------------------------


Pruning = Callable[..., bool]
SystemCondition = Callable[[PolynomialSystem], bool]


class Algorithm:
    def __init__(self, poly_system: PolynomialSystem,
                 system_conditions: Collection[SystemCondition] | None = None,
                 generation=default_generation,
                 scoring: Scoring = default_scoring,
                 pruning_funcs: Collection[Pruning] | None = None):
        self._system = poly_system
        self._generation = generation
        self._scoring = scoring
        self._pruning_funs = list(pruning_funcs) if pruning_funcs is not None else [lambda a, b, *_: False]
        self._sys_cond = list(system_conditions) if system_conditions else [lambda v: True]
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

    def quadratize(self) -> AlgorithmResult:
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
        return self.get_quadratizations(optimal_first.num_introduced_vars)

    @dump_results(log_enable, quad_systems_file)
    def get_quadratizations(self, depth: int) -> Set[PolynomialSystem]:
        return self.traverse_all(depth, lambda s: s.is_quadratized())

    def add_pruning(self, termination_criteria: Pruning) -> None:
        self._pruning_funs.append(termination_criteria)

    @property
    def generation_strategy(self):
        return self._generation

    @property
    def scoring(self):
        return self._scoring

    @generation_strategy.setter
    def generation_strategy(self, value):
        self._generation = value

    @logged(log_enable, log_file)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.generation_strategy. self.scoring)

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        pass


# ------------------------------------------------------------------------------

ALGORITHM_INTERRUPTED = False


def signal_handler(sig_num, frame):
    global ALGORITHM_INTERRUPTED
    print("The algorithm has been interrupted. Returning the current best.")
    ALGORITHM_INTERRUPTED = True


signal.signal(signal.SIGINT, signal_handler)


class BranchAndBound(Algorithm):

    def domination_upper_bound(self):
        system = self._system.copy()
        algo = BranchAndBound(system, self._sys_cond, self._generation, self._scoring,
                              [partial(pruning_by_domination, dominators=self.dominating_monomials),
                               *self._pruning_funs])
        res = algo.quadratize()
        upper_bound = res.num_introduced_vars
        if upper_bound != math.inf:
            self.preliminary_upper_bound = upper_bound + 1
        else:
            print("No upper bound was found")

    @timed(enabled=pb_enable)
    def quadratize(self) -> AlgorithmResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, self.preliminary_upper_bound)
        self._final_iter()
        self._save_results(opt_system)
        return AlgorithmResult(opt_system, nvars, traversed)

    @progress_bar(is_stop=False, enabled=pb_enable)
    def _bnb_step(self, part_res: PolynomialSystem, best_nvars) \
            -> Tuple[Union[int, float], Optional[PolynomialSystem], int]:
        self._nodes_traversed += 1
        # The order of these blocks is important: pruning rules assume that
        # the input partial result is not a quadratization
        if part_res.is_quadratized() and all(cond(part_res) for cond in self._sys_cond):
            return part_res.new_vars_count(), part_res, 1
        if any(map(lambda f: f(self, part_res, best_nvars), self._pruning_funs)) or ALGORITHM_INTERRUPTED:
            return math.inf, None, 1

        traversed_total = 1
        min_nvars, best_system = best_nvars, None
        for next_system in self.next_gen(part_res):
            nvars, opt_system, traversed = self._bnb_step(next_system, min_nvars)
            traversed_total += traversed
            if nvars < min_nvars:
                min_nvars = nvars
                best_system = opt_system
        return min_nvars, best_system, traversed_total

    @logged(log_enable, log_file)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.generation_strategy, self.scoring)

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @dump_results(enabled=log_enable, log_file=quad_systems_file)
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


def pruning_by_nodes_without_quadratization_found(algo: Algorithm, _: PolynomialSystem, *args, nodes_processed: int):
    best_nvars = args[0]
    if best_nvars < math.inf:
        return False
    elif algo._nodes_traversed >= nodes_processed:
        return True
    return False


def pruning_by_elapsed_time(algo: Algorithm, system: PolynomialSystem, *args, start_t, max_t):
    """
    Stops a search if 'max_t' was exceeded.

    :examples
        >>> from functools import partial
        >>> pruning = partial(pruning_by_elapsed_time, start_t=time(), max_t=100) # 100 seconds

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
            if part_res.laurent or all([x >= 0 for x in diff]):
                if diff in degree_one_monomials:
                    degree_one_monomials[diff] += 1
                else:
                    degree_one_monomials[diff] = 1
    degree_one_count = sorted(degree_one_monomials.values(), reverse=True)

    needed_new_vars = 0
    degree_two_monoms = len(part_res.nonsquares)
    while True:
        if degree_two_monoms <= (needed_new_vars * (needed_new_vars + 1)) // 2:
            break
        if needed_new_vars < len(degree_one_count):
            degree_two_monoms -= degree_one_count[needed_new_vars]
        needed_new_vars += 1

    lower_bound = part_res.new_vars_count() + needed_new_vars

    if lower_bound >= best_nvars:
        return True
    return False


def pruning_by_declining_variables(a: Algorithm, part_res: PolynomialSystem, *args, excl_vars: List[Tuple]):
    return not without_variables(part_res, excl_vars)


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
            if part_res.laurent or all([x >= 0 for x in diff]):
                if diff in degree_one_monomials:
                    degree_one_monomials[diff] += 1
                else:
                    degree_one_monomials[diff] = 1
    degree_one_count = sorted(degree_one_monomials.values(), reverse=True)

    needed_new_vars = 0
    while True:
        if needed_new_vars >= len(MAX_C4_FREE_EDGES):
            break
        if no_C4_edges <= MAX_C4_FREE_EDGES[needed_new_vars][min(max_loops, needed_new_vars)]:
            break
        if needed_new_vars < len(degree_one_count):
            no_C4_edges -= degree_one_count[needed_new_vars]
        needed_new_vars += 1

    lower_bound = part_res.new_vars_count() + needed_new_vars

    if lower_bound >= best_nvars:
        return True

    return False


def pruning_by_domination(a: BranchAndBound, part_res: PolynomialSystem, *args, dominators):
    if not part_res.introduced_vars:
        return False

    return not all([dominated(m, dominators) for m in part_res.introduced_vars])


default_pruning_rules = [
    pruning_by_quadratic_upper_bound,
    pruning_by_squarefree_graphs
]


# ------------------------------------------------------------------------------------------------

def without_variables(part_res: PolynomialSystem, excl_vars: List[Tuple]):
    excl_indices = [np.argmax(v) for v in excl_vars]
    if any([any([v[i] != 0 for i in excl_indices]) for v in part_res.introduced_vars]):
        return False
    return True


def with_higher_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) > system.original_degree, system.vars))


def with_le_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) <= system.original_degree, system.vars))
