import math
import itertools
import configparser
import pickle
from sympy import *
from collections import deque
from typing import Callable, List, Optional, Generator, Set
from functools import partial, reduce
from operator import mul
from random import randrange
from .heuristics import *  # replace with .heuristics if you want pip install
from .util import *  # replace with .util if you want pip install

from memory_profiler import profile

config = configparser.ConfigParser({'logging_enable': False, 'logging_file': '../log/log.csv'})
config.read("../config.ini")
log_enable = eval(config.get('DEFAULT', 'logging_enable'))  # Security Error here, but does not matter I believe
log_file = config.get('DEFAULT', 'logging_file')


# ------------------------------------------------------------------------------

class PolynomialSystem:
    def __init__(self, polynomials: List[sp.Poly]):
        """
        polynomials - right-hand sides of the ODE system listed in the same order as
                      the variables in the polynomial ring
        """
        self.dim = len(polynomials[0].ring.gens)
        self.gen_syms = list(map(lambda g: sp.Symbol(str(g)), polynomials[0].ring.gens))

        # put not monomials but differents in the exponents between rhs and lhs
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

    def add_var(self, v):
        for i in range(self.dim):
            if v[i] > 0:
                for m in self.rhs[i]:
                    self.nonsquares.add(tuple([v[j] + m[j] for j in range(self.dim)]))

        self.vars.add(v)
        for u in self.vars:
            self.squares.add(tuple([u[i] + v[i] for i in range(self.dim)]))
        self.nonsquares = set(filter(lambda s: s not in self.squares, self.nonsquares))

    def is_quadratized(self):
        return not self.nonsquares

    def get_smallest_nonsquare(self):
        return min([(sum(m), m) for m in self.nonsquares])[1]

    def next_generation(self, heuristics=default_score):
        if len(self.nonsquares) == 0:
            return list()
        new_gen = []
        for d in get_decompositions(self.get_smallest_nonsquare()):
            c = pickle.loads(pickle.dumps(self, -1))
            for v in d:
                c.add_var(v)
            new_gen.append(c)

        return sorted(new_gen, key=heuristics)

    def new_vars_count(self):
        return len(self.vars) - self.dim - 1

    def apply_quadratizaton(self):
        # TODO
        new_variables = list(filter(lambda m: monomial_deg(m) >= 2, self.vars))

    def to_sympy(self, gens):
        # TODO: Does not work correctly with PolyElement
        return [mlist_to_poly(mlist, gens) for mlist in self.rhs.values()]

    def __repr__(self):
        return f"{self._introduced_variables_str()}"

    def _introduced_variables_str(self):
        return list(map(lambda v: latex(monomial_to_poly(Monomial(v, self.gen_syms)).as_expr()),
                        self._remove_free_variables(self.vars)))

    def _remove_free_variables(self, vars: Collection[tuple]):
        return tuple(filter(lambda v: monomial_deg(v) >= 2, vars))


# ------------------------------------------------------------------------------

class QuadratizationResult:
    def __init__(self,
                 system: Optional[PolynomialSystem],
                 introduced_vars: int,
                 nodes_traversed: int):
        self.system = system
        self.introduced_vars = introduced_vars
        self.nodes_traversed = nodes_traversed

    def sympy_str(self, variables):
        pass  # TODO: apply quadratization to system

    def __repr__(self):
        if self.system is None:
            return "No quadratization found under the given condition\n" + \
                   f"Nodes traversed: {self.nodes_traversed}"
        return f"Number of introduced variables: {self.introduced_vars}\n" + \
               f"Introduced variables: {self.system._introduced_variables_str()}\n" + \
               f"Nodes traversed: {self.nodes_traversed}"


# ------------------------------------------------------------------------------

EarlyTermination = Callable[..., bool]


class Algorithm:
    def __init__(self,
                 poly_system: PolynomialSystem,
                 heuristics: Heuristics = default_score,
                 early_termination: Collection[EarlyTermination] = None):
        self._system = poly_system
        self._heuristics = heuristics
        self._early_termination_funs = list(early_termination) if early_termination is not None else [
            lambda a, b, *_: False]
        self._nodes_traversed = 0

    def quadratize(self) -> QuadratizationResult:
        pass

    def traverse_all(self, to_depth: int, pred: Callable[[PolynomialSystem], bool]):
        res = set()
        self._dls(self._system, to_depth, pred, res)
        self._final_iter()
        return res

    @progress_bar(is_stop=False)
    def _dls(self, part_res: PolynomialSystem, to_depth: int, pred: Callable[[PolynomialSystem], bool], res: set):
        if part_res.new_vars_count() > to_depth:
            return

        if pred(part_res):
            res.add(part_res)
        else:
            for next_system in part_res.next_generation():
                self._dls(next_system, to_depth, pred, res)
        return

    def get_optimal_quadratizations(self) -> Set[PolynomialSystem]:
        optimal_first = self.quadratize()
        print(optimal_first)
        return self.get_quadratizations(optimal_first.introduced_vars)

    def get_quadratizations(self, depth: int) -> Set[PolynomialSystem]:
        return self.traverse_all(depth, lambda s: s.is_quadratized())

    def attach_early_termimation(self, termination_criteria: EarlyTermination) -> None:
        self._early_termination_funs.append(termination_criteria)

    @property
    def heuristics(self):
        return self._heuristics

    @heuristics.setter
    def heuristics(self, value):
        self._heuristics = value

    @progress_bar(is_stop=True)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        pass


# ------------------------------------------------------------------------------

class BranchAndBound(Algorithm):
    def __init__(self, poly_system: PolynomialSystem,
                 heuristics: Heuristics = default_score,
                 early_termination: Union[EarlyTermination, Collection[EarlyTermination]] = None):
        super().__init__(poly_system, heuristics, early_termination)

    @timed
    def quadratize(self, cond: Callable[[PolynomialSystem], bool] = lambda _: True) -> QuadratizationResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, math.inf, cond)
        self._final_iter()
        return QuadratizationResult(opt_system, nvars, traversed)

    @progress_bar(is_stop=False)
    def _bnb_step(self, part_res: PolynomialSystem, best_nvars, cond) \
            -> Tuple[Union[int, float], Optional[PolynomialSystem], int]:
        self._nodes_traversed += 1
        if part_res.is_quadratized() and cond(part_res):
            return part_res.new_vars_count(), part_res, 1
        if any(map(lambda f: f(self, part_res, best_nvars), self._early_termination_funs)):
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
        return part_res.next_generation(self.heuristics)

    @progress_bar(is_stop=True)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0


# ------------------------------------------------------------------------------

class ID_DLS(Algorithm):
    def __init__(self, poly_system: PolynomialSystem,
                 start_upper_bound: int,
                 upper_bound: int,
                 heuristics: Heuristics = default_score,
                 early_termination: Union[EarlyTermination, Collection[EarlyTermination]] = None):
        super().__init__(poly_system, heuristics, early_termination)
        self.upper_bound = upper_bound
        self.start_upper_bound = start_upper_bound

    @timed
    def quadratize(self) -> QuadratizationResult:
        stack = deque()
        high_depth_stack = deque()

        curr_depth = 0
        curr_max_depth = self.start_upper_bound
        stack.append((self._system, curr_depth))

        while True:
            self._iter()
            self._nodes_traversed += 1
            if len(stack) == 0:
                if len(high_depth_stack) == 0:
                    raise RuntimeError("Limit depth passed. No quadratic system is found.")
                stack = high_depth_stack
                high_depth_stack = deque()
                curr_max_depth += int(math.ceil(math.log(curr_depth + 1)))

            system, curr_depth = stack.pop()
            if (any(map(lambda f: f(self, system), self._early_termination_funs))):
                return QuadratizationResult(None, -1, self._nodes_traversed)
            if system.is_quadratized():
                traversed = self._nodes_traversed
                self._final_iter()
                return QuadratizationResult(system, curr_depth, traversed)

            if curr_depth > self.upper_bound:
                continue

            for next_system in system.next_generation(self.heuristics):
                if curr_depth < curr_max_depth:
                    stack.append((next_system, curr_depth + 1))
                else:
                    high_depth_stack.append((next_system, curr_depth + 1))

    @progress_bar(is_stop=True)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @progress_bar(is_stop=False)
    def _iter(self):
        pass


# ------------------------------------------------------------------------------

def termination_by_nodes_processed(algo: Algorithm, _: PolynomialSystem, *args, nodes_processed: int):
    if algo._nodes_traversed >= nodes_processed:
        return True
    return False


def termination_by_vars_number(_: Algorithm, system: PolynomialSystem, *args, nvars: int):
    if system.new_vars_count() >= nvars:
        return True
    return False


def termination_by_best_nvars(a: Algorithm, part_res: PolynomialSystem, *args):
    best_nvars, *_ = args
    if part_res.new_vars_count() >= best_nvars - 1:
        return True
    return False


def termination_by_square_bound(a: Algorithm, part_res: PolynomialSystem, *args):
    best_nvars, *_ = args
    total_monoms = len(part_res.squares) + len(part_res.nonsquares)
    lower_bound = int(math.ceil((math.sqrt(1. + 8. * total_monoms) - 1.) / 2.))
    if lower_bound >= best_nvars - 1:  # TODO: check correctness here
        return True
    return False


def termination_by_C4_bound(a: Algorithm, part_res: PolynomialSystem, *args):
    best_nvars, *_ = args
    no_C4_monoms = set()
    sums_of_monoms = set()
    for m in sorted(part_res.squares.union(part_res.nonsquares), key=sum, reverse=True):
        new_sums = set()
        to_add = True
        for mm in no_C4_monoms:
            s = tuple([m[i] + mm[i] for i in range(len(m))])
            if s in sums_of_monoms:
                to_add = False
                break
            new_sums.add(s)
        if to_add:
            sums_of_monoms = sums_of_monoms.union(new_sums)
            no_C4_monoms.add(m)

    m = len(no_C4_monoms)
    lb = 1
    while 4 * m ** 2 - 2 * lb * m - lb ** 2 * (lb + 1) > 0:
        lb += 1
    if lb >= best_nvars - 1:  # TODO: check correctness here
        return True
    return False


def with_higher_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) > system.original_degree, system.vars))


def with_le_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) <= system.original_degree, system.vars))
