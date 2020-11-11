import math
from sympy import *
from collections import deque
from typing import Callable, Union, Optional
from heuristics import *
from structures import PolynomialSystem
from util import *


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


if __name__ == "__main__":
    R, x = ring(["x"], QQ)
    poly_system = PolynomialSystem([x ** 6])
    algo = BranchAndBound(poly_system, 3)
    res = algo.quadratize()
    print(res)
