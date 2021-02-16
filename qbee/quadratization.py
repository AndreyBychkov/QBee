import math
import itertools
import hashlib
import configparser
import pickle
import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from sympy import *
from sympy.polys.orderings import monomial_key
from collections import deque
from typing import Callable, List, Optional, Generator, Set
from functools import partial, reduce
from operator import mul
from random import randrange
from .heuristics import *  # replace with .heuristics if you want pip install
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

    @property
    def quadratization(self):
        return sorted(map(lambda v: monomial_to_poly(Monomial(v, self.gen_syms)).as_expr(),
                          self.introduced_vars), key=monomial_key('grlex', self.gen_syms))

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
        return min([(sum(m), m) for m in self.nonsquares])[1]

    def next_generation(self, heuristics=default_score):
        if len(self.nonsquares) == 0:
            return list()
        new_gen = []
        for d in get_decompositions(self.get_smallest_nonsquare()):
            c = pickle.loads(pickle.dumps(self, -1))  # inline self.copy for speedup
            for v in d:
                c.add_var(v)
            new_gen.append(c)

        return sorted(new_gen, key=heuristics)

    def new_vars_count(self):
        return len(self.vars) - self.dim - 1

    def to_sympy(self, gens):
        # TODO: Does not work correctly with PolyElement
        return [mlist_to_poly(mlist, gens) for mlist in self.rhs.values()]

    def apply_quadratization(self):
        """TODO: Not fully done yet"""

        def add_introduced_equations():
            rhs = self.rhs.copy()
            n0 = len(rhs)

            pass

        qvars = self._name_vars()
        rhs = [list(v) for k, v in sorted(self.rhs.items(), key=lambda kv: kv[0])]

        def quad_var_change(monom):
            vs = self.vars
            for v1 in vs:
                for v2 in vs:
                    if np.array_equal(np.array(v1) + np.array(v2), monom):
                        return (v1, v2)
            raise Exception("Monomial can not be represented as quadratic, but quadratization is found."
                            " Contact developer team")

        def quad_monom_to_latex(monom):
            def to_symbol(m):
                if any(m):
                    return sp.Symbol(qvars[m])
                return sp.Integer(1)

            return latex(reduce(lambda a, b: a * b, [to_symbol(m) for m in sorted(monom)]))

        for i, eq in enumerate(rhs):
            for j, m in enumerate(eq):
                m = list(m)
                m[i] += 1
                m = quad_var_change(m)
                eq[j] = quad_monom_to_latex(m)
        return reduce(lambda a, b: a + ' \n' + b, [reduce(lambda x, y: x + ' + ' + y, eq) for eq in rhs])

    def _name_vars(self) -> dict:
        return {**self._name_orig_vars(), **self._name_introduced_vars(), **self._name_constant_var()}

    def _name_orig_vars(self):
        def make_tuple(i: int):
            t = [0, ] * len(self.gen_syms)
            t[i] = 1
            return tuple(t)

        return dict([(make_tuple(i), latex(s)) for i, s in enumerate(self.gen_syms)])

    def _name_constant_var(self):
        return {(0,) * len(self.gen_syms): sp.Integer(1)}

    def _name_introduced_vars(self):
        base = 'w'
        return dict([(v, f"{base}{i + 1}") for i, v in enumerate(self.introduced_vars)])

    def __str__(self):
        return f"{self._introduced_variables_str()}"

    def __repr__(self):
        return f"{self._introduced_variables_str()}"

    def _introduced_variables_str(self):
        """Legacy compatibility str representation. TODO: refactor it"""
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

EarlyTermination = Callable[..., bool]


class Algorithm:
    def __init__(self,
                 poly_system: PolynomialSystem,
                 heuristics: Heuristics = default_score,
                 early_termination: Collection[EarlyTermination] = None,
                 use_weak_hull=True):
        self._system = poly_system
        self._heuristics = heuristics
        self._early_termination_funs = list(early_termination) if early_termination is not None else [
            lambda a, b, *_: False]
        self._nodes_traversed = 0

        self.hull: Optional[ConvexHull] = None
        if len(list(poly_system.vars)[0]) > 1:
            points = list()
            for i, eq in poly_system.rhs.items():
                for m in eq:
                    new_m = list(m)
                    new_m[i] = new_m[i] if new_m[i] > 0 else 0
                    points.append(tuple(new_m))
            self.hull = ConvexHull(points + list(poly_system.vars))

        self.weak_hull: Optional[ConvexHull] = None
        if self.hull and use_weak_hull:
            points = [(0, ) * poly_system.dim, ]
            max_order = max(list(map(sum, self.hull.points[self.hull.vertices].astype(int).tolist())))
            for v in range(poly_system.dim):
                p = [0] * poly_system.dim
                p[v] = max_order
                points.append(tuple(p))
            self.weak_hull = ConvexHull(points)
            self.attach_early_termination(partial(termination_by_newton_polyhedron, hull=Delaunay(self.weak_hull.points)))
            print(f"Weak convex hull is is set to order = {max_order}")


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

    def attach_early_termination(self, termination_criteria: EarlyTermination) -> None:
        self._early_termination_funs.append(termination_criteria)

    @property
    def heuristics(self):
        return self._heuristics

    @heuristics.setter
    def heuristics(self, value):
        self._heuristics = value

    @logged(log_enable, log_file)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.heuristics)

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        pass


# ------------------------------------------------------------------------------

class BranchAndBound(Algorithm):

    def newton_polyhedral_vertices_upper_bound(self):
        system = self._system.copy()
        hull_vars = list(map(tuple, self.hull.points[self.hull.vertices].astype(int)))
        hull_vars = list(filter(lambda v: v not in system.vars, hull_vars))
        for v in hull_vars:
            system.add_var(v)
        algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
        quad_res = algo.quadratize()
        upper_bound = quad_res.introduced_vars
        if upper_bound != math.inf:
            self.attach_early_termination(partial(termination_by_vars_number, nvars=upper_bound))
            print(f"Upper bound for quadratization is {upper_bound}")
        else:
            print(f"Upper bound is not found")

    def inside_newton_polyhedral_upper_bound(self):
        system = self._system.copy()
        algo = BranchAndBound(system, aeqd_score,
                              [termination_by_best_nvars,
                               partial(termination_by_newton_polyhedron, hull=Delaunay(self.hull.points))])
        res = algo.quadratize()
        upper_bound = res.introduced_vars
        if upper_bound != math.inf:
            self.attach_early_termination(partial(termination_by_vars_number, nvars=upper_bound))
            print(f"Upper bound for quadratization is {upper_bound}")
        else:
            print(f"Upper bound is not found")

    @timed
    def quadratize(self, cond: Callable[[PolynomialSystem], bool] = lambda _: True) -> QuadratizationResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, math.inf, cond)
        self._final_iter()
        self._save_results(opt_system)
        return QuadratizationResult(opt_system, nvars, traversed)

    @progress_bar(is_stop=False, enabled=pb_enable)
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

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @dump_results(enabled=True, log_file=quad_systems_file)
    def _save_results(self, opt_system):
        return [opt_system, ]


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

    @progress_bar(is_stop=True, enabled=pb_enable)
    @logged(log_enable, log_file, is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @progress_bar(is_stop=False, enabled=pb_enable)
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


def termination_by_newton_polyhedron(a: BranchAndBound, part_res: PolynomialSystem, *args, hull: Delaunay):
    if not part_res.introduced_vars:
        return False

    return not all(hull.find_simplex(part_res.introduced_vars) >= 0)


def with_higher_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) > system.original_degree, system.vars))


def with_le_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) <= system.original_degree, system.vars))
