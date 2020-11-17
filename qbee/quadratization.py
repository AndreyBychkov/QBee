import math
import copy
from sympy import *
from collections import deque
from typing import Callable, List, Optional, Generator, Set
from functools import partial, reduce
from itertools import product
from operator import mul
from random import randrange
from heuristics import *  # replace with .heuristics if you want pip install
from util import *  # replace with .util if you want pip install


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
        return [mlist_to_poly(mlist, gens) for mlist in self.rhs.values()]

    def __repr__(self):
        return f"{self._intoroduced_variables_str()}"

    def _intoroduced_variables_str(self):
        return list(map(lambda v: latex(monomial_to_poly(Monomial(v, self.gen_syms)).as_expr()),
                        self._remove_free_variables(self.vars)))

    def _remove_free_variables(self, vars: Collection[tuple]):
        return tuple(filter(lambda v: monomial_deg(v) >= 2, vars))


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
               f"Introduced variables: {self.system._intoroduced_variables_str()}\n" + \
               f"Nodes traversed: {self.nodes_traversed}"


EarlyTermination = Callable[..., bool]


class Algorithm:
    def __init__(self,
                 poly_system: PolynomialSystem,
                 heuristics: Heuristics = default_score,
                 early_termination: Collection[EarlyTermination] = None):
        self._system = poly_system
        self._heuristics = heuristics
        self._early_termination_funs = list(early_termination) if early_termination is not None else [
            lambda a, b: False]
        self._nodes_traversed = 0

    def quadratize(self) -> QuadratizationResult:
        pass

    def traverse_all(self, to_depth: int, pred: Callable[[PolynomialSystem], bool]):
        res = set()
        self._dls(self._system, to_depth, pred, res)
        self._final_iter()
        return res

    @logged(is_stop=False)
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
        return self.traverse_all(optimal_first.introduced_vars, lambda s: s.is_quadratized())

    def attach_early_termimation(self, termination_criteria: EarlyTermination) -> None:
        self._early_termination_funs.append(termination_criteria)

    @property
    def heuristics(self):
        return self._heuristics

    @heuristics.setter
    def heuristics(self, value):
        self._heuristics = value

    @logged(is_stop=True)
    def _final_iter(self):
        pass


class BranchAndBound(Algorithm):
    def __init__(self, poly_system: PolynomialSystem,
                 heuristics: Heuristics = default_score,
                 early_termination: Union[EarlyTermination, Collection[EarlyTermination]] = None):
        super().__init__(poly_system, heuristics, early_termination)

    @timed
    def quadratize(self) -> QuadratizationResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, math.inf)
        self._final_iter()
        return QuadratizationResult(opt_system, nvars, traversed)

    @logged(is_stop=False)
    def _bnb_step(self, part_res: PolynomialSystem, best_nvars) \
            -> Tuple[Union[int, float], Optional[PolynomialSystem], int]:
        self._nodes_traversed += 1
        if part_res.is_quadratized():
            return part_res.new_vars_count(), part_res, 1
        if any(map(lambda f: f(self, part_res, best_nvars), self._early_termination_funs)):
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

    @logged(is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0


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

    @logged(is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @logged(is_stop=False)
    def _iter(self):
        pass


def termination_by_nodes_processed(algo: Algorithm, _: PolynomialSystem, *args, nodes_processed: int):
    if algo._nodes_traversed >= nodes_processed:
        return True
    return False


def termination_by_vars_number(_: Algorithm, system: PolynomialSystem, *args, nvars: int):
    if len(system.vars) >= nvars:
        return True
    return False


def termination_by_best_nvars(a: Algorithm, part_res: PolynomialSystem, *args):
    best_nvars, *_ = args
    if part_res.new_vars_count() >= best_nvars - 1:
        return True
    return False


class SimpleGenerator:
    def __init__(self, variables_num, system_degree):
        self.nvars = variables_num
        self.deg = system_degree
        _, *self.vars = ring([f"x{i}" for i in range(self.nvars)], QQ)

    def generate(self):
        return [self._gen_eq() for _ in range(self.nvars)]

    def _gen_eq(self):
        return reduce(add, [self._gen_eq_deg(d) for d in range(1, self.deg + 1)])

    def _gen_eq_deg(self, deg: int):
        return reduce(add, [randrange(-10, 10) * reduce(mul, p) for p in product(self.vars, repeat=deg)]) + randrange(
            -10, 10)


def run_with_gen(N):
    print(N)
    gen = SimpleGenerator(1, N)
    system = gen.generate()
    poly_system = PolynomialSystem(system)
    algo = BranchAndBound(poly_system, 13, heuristics=default_score)
    res = algo.quadratize()
    print(res)


if __name__ == "__main__":
    R, x = ring(["x", ], QQ)
    poly_system = PolynomialSystem([x ** 10 + x ** 9 + x ** 8 + x ** 5 + x + 1])
    algo = BranchAndBound(poly_system, heuristics=aeqd_score,
                          early_termination=[termination_by_best_nvars])
    for res in algo.get_optimal_quadratizations():
        print(res)
