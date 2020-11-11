import copy
import math
import time
import numpy as np

from sympy import *
from sympy.polys.monomials import monomial_deg, monomial_divides
from typing import Callable, Union, Collection, Optional, List, Tuple
from functools import reduce
from operator import add
from util import monomial_to_poly
from tqdm import tqdm
from collections import deque


def timed(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        res = func(*args, **kwargs)
        end_time = time.time()
        print()
        print(f"Elapsed time: {np.round(end_time - start_time, 3)}s.")
        return res

    return wrapper


def get_decompositions(monomial):
    if len(monomial) == 0:
        return {(tuple(), tuple())}
    result = set()
    prev_result = get_decompositions(tuple(monomial[:-1]))
    for r in prev_result:
        for i in range(monomial[-1] + 1):
            a, b = tuple(list(r[0]) + [i]), tuple(list(r[1]) + [monomial[-1] - i])
            result.add((min(a, b), max(a, b)))
    return result


def empty_score(system) -> int:
    return 1


def default_score(system) -> int:
    total_nonsquare = sum([monomial_deg(m) for m in system.nonsquares])
    return total_nonsquare + system.dim * len(system.vars)


def aeqd_score(system) -> int:
    eq_degs = list(map(lambda mlist: max(map(monomial_deg, mlist)), system.rhs.values()))
    aeqds = list(map(_compute_aeqd, system.nonsquares, [eq_degs, ] * len(system.nonsquares)))
    return sum(aeqds)


def smd_score(system) -> int:
    mlist = system.nonsquares
    return sum(map(lambda s: _compute_smd(s, mlist), mlist))


def _compute_aeqd(sub: Tuple[int], eq_degs):
    mon_degs = map(lambda deg: deg + monomial_deg(sub) - 1, eq_degs)
    quad_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, mon_degs))
    return sum(quad_discrepancies)


def _compute_smd(sub, mlist: list):
    return (monomial_deg(sub) - 1) * len(list(filter(lambda m: monomial_divides(sub, m), mlist)))


def _mlist_to_poly(mlist: Collection[Monomial], gens) -> Poly:
    return Poly(reduce(add, mlist), gens)


class PolynomialSystem:
    def __init__(self, polynomials: List[Poly]):
        """
        polynomials - right-hand sides of the ODE system listed in the same order as 
                      the variables in the polynomial ring
        """
        self.dim = len(polynomials[0].ring.gens)
        self.gen_syms = list(map(lambda g: Symbol(str(g)), polynomials[0].ring.gens))

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
        new_gen = []
        for d in get_decompositions(self.get_smallest_nonsquare()):
            c = copy.deepcopy(self)
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
        return [_mlist_to_poly(mlist, gens) for mlist in self.rhs.values()]


Heuristics = Callable[[PolynomialSystem], int]
TerminationCriteria = Callable[[PolynomialSystem], bool]


class QuadratizationResult:
    def __init__(self,
                 system: PolynomialSystem,
                 introduced_vars: int,
                 nodes_traversed: int):
        self.system = system
        self.introduced_vars = introduced_vars
        self.nodes_traversed = nodes_traversed

    def sympy_str(self, variables):
        pass

    def __repr__(self):
        return f"Number of introduced variables: {self.introduced_vars}\n" + \
               f"Introduced variables: {self._intoroduced_variables_str()}\n" + \
               f"Nodes traversed: {self.nodes_traversed}"

    def _intoroduced_variables_str(self):
        return list(map(lambda v: latex(monomial_to_poly(Monomial(v, self.system.gen_syms)).as_expr()),
                        self._remove_free_variables(self.system.vars)))

    def _remove_free_variables(self, vars: Collection[tuple]):
        return tuple(filter(lambda v: monomial_deg(v) >= 2, vars))

    def _var_to_symbolic(self, v: tuple):
        return ''.join([f"{g}^{p}" for g, p in zip(self.system.gen_syms, v)])


class Algorithm:
    def __init__(self,
                 poly_system: PolynomialSystem,
                 heuristics: Heuristics = default_score,
                 termination_criteria: Union[TerminationCriteria, Collection[TerminationCriteria]] = None):
        self._system = poly_system
        self._heuristics = heuristics
        self._early_termination_funs = list(termination_criteria) if termination_criteria is not None else [
            lambda _: False, ]

    def quadratize(self) -> QuadratizationResult:
        pass

    def attach_early_termimation(self, termination_criteria: Callable[[PolynomialSystem], bool]) -> None:
        self._early_termination_funs.append(termination_criteria)

    @property
    def heuristics(self):
        return self._heuristics

    @heuristics.setter
    def heuristics(self, value):
        self._heuristics = value


class BranchAndBound(Algorithm):
    def __init__(self, poly_system: PolynomialSystem,
                 upper_bound: int,
                 heuristics: Heuristics = default_score,
                 termination_criteria: Union[TerminationCriteria, Collection[TerminationCriteria]] = None):
        super().__init__(poly_system, heuristics, termination_criteria)
        self.upper_bound = upper_bound

    def quadratize(self) -> QuadratizationResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, self.upper_bound)
        return QuadratizationResult(opt_system, nvars, traversed)

    def _bnb_step(self, part_res: PolynomialSystem, best_nvars) -> Tuple[
        Union[int, float], Optional[PolynomialSystem], int]:
        if part_res.is_quadratized():
            return part_res.new_vars_count(), part_res, 1
        if part_res.new_vars_count() >= best_nvars - 1:
            return math.inf, None, 1

        traversed_total = 1
        min_nvars, best_system = best_nvars, None
        for next_system in part_res.next_generation(self.heuristics):
            nvars, opt_system, traversed = self._bnb_step(next_system, min_nvars)
            traversed_total += traversed
            if nvars < min_nvars:
                min_nvars = nvars
                best_system = opt_system
        return min_nvars, best_system, traversed_total


class ID_DLS(Algorithm):
    def __init__(self, poly_system: PolynomialSystem,
                 start_upper_bound: int,
                 upper_bound: int,
                 heuristics: Heuristics = default_score,
                 termination_criteria: Union[TerminationCriteria, Collection[TerminationCriteria]] = None):
        super().__init__(poly_system, heuristics, termination_criteria)
        self.upper_bound = upper_bound
        self.start_upper_bound = start_upper_bound

    def quadratize(self) -> QuadratizationResult:
        stack = deque()
        high_depth_stack = deque()

        curr_depth = 0
        curr_max_depth = self.start_upper_bound
        stack.append((self._system, curr_depth))
        nodes_traversed = 1

        while True:
            if len(stack) == 0:
                if len(high_depth_stack) == 0:
                    raise RuntimeError("Limit depth passed. No quadratic system is found.")
                stack = high_depth_stack
                high_depth_stack = deque()
                curr_max_depth += int(math.ceil(math.log(curr_depth + 1)))

            system, curr_depth = stack.pop()
            nodes_traversed += 1
            if system.is_quadratized():
                return QuadratizationResult(system, curr_depth, curr_depth)

            if curr_depth > self.upper_bound:
                continue

            for next_system in system.next_generation(self.heuristics):
                if curr_depth < curr_max_depth:
                    stack.append((next_system, curr_depth + 1))
                else:
                    high_depth_stack.append((next_system, curr_depth + 1))


class BestFirst(Algorithm):
    def __init__(self, poly_system: PolynomialSystem,
                 upper_bound: int,
                 heuristics: Heuristics = default_score,
                 termination_criteria: Union[TerminationCriteria, Collection[TerminationCriteria]] = None):
        super().__init__(poly_system, heuristics, termination_criteria)
        self.upper_bound = upper_bound

    def quadratize(self) -> QuadratizationResult:
        pass


@timed
def simple_test():
    R, x, y = ring(['x', 'y'], QQ)
    poly_system = PolynomialSystem([x ** 2 + y, 3 * x * y ** 3 - y])
    algo = BranchAndBound(poly_system, 10, heuristics=empty_score)
    res = algo.quadratize()
    print(res)


@timed
def poly_test():
    R, x, y = ring(['x', 'y'], QQ)
    poly_system = PolynomialSystem([(x + 1) ** 8, y])
    algo = BranchAndBound(poly_system, 20, heuristics=empty_score)
    res = algo.quadratize()
    print(res)


@timed
def long_poly_test():
    R, x, y = ring(['x', 'y'], QQ)
    poly_system = PolynomialSystem([(x + 1) ** 15, y])
    algo = BranchAndBound(poly_system, 8, heuristics=empty_score)
    res = algo.quadratize()
    print(res)


@timed
def xSigmoid():
    R, x, y, z = ring(["x", "y", "z"], QQ)
    poly_system = PolynomialSystem([x * z, x * y * z, x * y * z ** 3])
    algo = BranchAndBound(poly_system, 5, heuristics=default_score)
    res = algo.quadratize()
    print(res)


@timed
def RabinovichFabricant():
    R, x, y, z = ring(["x", "y", "z"], QQ)
    poly_system = PolynomialSystem([y * (z - 1 - x ** 2) + x, x * (3 * z + 1 - x ** 2) + y, -2 * z * (2 + x * y)])
    algo = BranchAndBound(poly_system, 20, heuristics=smd_score)
    res = algo.quadratize()
    print(res)


if __name__ == "__main__":
    R, x = ring(["x"], QQ)
    poly_system = PolynomialSystem([x ** 6])
    algo = BranchAndBound(poly_system, 3)
    res = algo.quadratize()
    print(res)
