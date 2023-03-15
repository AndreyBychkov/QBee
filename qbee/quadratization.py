from __future__ import annotations

import copy
import math
import signal
import pickle
import numpy as np
import sympy as sp
from sympy.polys.rings import PolyElement
from typing import Set, Collection
from ordered_set import OrderedSet
from functools import partial
from time import time
from .selection import *
from .util import progress_bar, dump_results, dominated, monom2PolyElem, apply_quadratization, monom2str, derivatives, \
    logged
from .polynomialization import EquationSystem, polynomialize, eq_list_to_eq_system
from .printer import str_qbee


def polynomialize_and_quadratize(system: EquationSystem | list[(sp.Symbol, sp.Expr)],
                                 input_free=False, input_der_orders=None,
                                 conditions: Collection["SystemCondition"] = (),
                                 polynomialization_upper_bound=10, calc_upper_bound=True,
                                 generation_strategy=default_generation, scoring: Scoring = default_scoring,
                                 pruning_functions: Collection["Pruning"] | None = None,
                                 new_vars_name="w_", start_new_vars_with=0) -> QuadratizationResult | None:
    """
    Polynomialize and then quadratize a system of ODEs with the continuous right-hand side.

    :param system: system of equations in form [(X, f(X)), ...] where the left-hand side is derivatives of X.
    :param input_free: if True the function will not introduce derivatives of input functions.
    :param input_der_orders: a mapping of input variables to maximum order of their derivatives.
     For example, {T: 2} => T' and T'' are introduced.
    :param conditions: a list of predicates PolynomialSystem -> bool the quadratized systems must comply with.
     For example, use `partial(without_variables, [x1, x2])` to stop the algorithm from using x1 and x2 in new varaibles.
    :param polynomialization_upper_bound: how many new variables the polynomialization algorithm can introduce.
    :param calc_upper_bound: if True a non-optimal quadratization will be quickly found and used as upper bound.
     Disable if you already know an upper bound and add it in pruning_functions via `pruning_by_vars_number`.
    :param generation_strategy: a function for proposing new variables during quadratization evaluation.
    :param pruning_functions: a predicate indicating should the algorithm drop a current branch.
     Unlike `conditions` parameter, prunings apply to each intermediate result, not just to quadratizations.
     Check functions starting with `pruning_by` to see examples.
    :param new_vars_name: base name for new variables. For example, new_var_name='z' => z0, z1, z2, ...
    :param start_new_vars_with: initial index for new variables. Example: start_new_vars_with=3 => w_3, w_4, ...
    :return: a container of a quadratized system and new variables introduced or None if there is nothing found

    Example:
        >>> from qbee import *
        >>> from sympy import exp
        >>> x, y, u = functions("x, y, u")
        >>> p = parameters("p")
        >>> polynomialize_and_quadratize([(x, y / (1 + exp(-p * x))), (y, x * exp(y) + u)], input_free=True).print()
        Introduced variables:
        w_0 = 1/(1 + exp(-p*x))
        w_1 = exp(y)
        w_2 = exp(-p*x)
        w_3 = w_1*x
        w_4 = w_0*y
        w_5 = w_0**2*w_2*y
        w_6 = w_0**2*w_2
         ‎
        x' = w_4
        y' = u + w_3
        w_0' = p*w_4*w_6
        w_1' = u*w_1 + w_1*w_3
        w_2' = -p*w_2*w_4
        w_3' = u*w_3 + w_1*w_4 + w_3**2
        w_4' = p*w_4*w_5 + u*w_0 + w_0*w_3
        w_5' = -p*w_4*w_5 + 2*p*w_5**2 + u*w_6 + w_3*w_6
        w_6' = -p*w_4*w_6 + 2*p*w_5*w_6

        >>> polynomialize_and_quadratize([(x, y**3), (y, x**3)], new_vars_name="c_", start_new_vars_with=1).print()
        Introduced variables:
        c_1 = y**2
        c_2 = x**2
        c_3 = x*y
         ‎
        x' = c_1*y
        y' = c_2*x
        c_1' = 2*c_2*c_3
        c_2' = 2*c_1*c_3
        c_3' = c_1**2 + c_2**2

        >>> upper_bound = partial(pruning_by_vars_number, nvars=10)
        >>> res = polynomialize_and_quadratize([(x, x**2 * u)], input_free=True, pruning_functions=[upper_bound, *default_pruning_rules])
        >>> print(res is None)
        True
    """
    if input_der_orders is None:
        input_der_orders = dict()
    poly_system = polynomialize(system, polynomialization_upper_bound,
                                new_vars_name=new_vars_name, start_new_vars_with=start_new_vars_with)
    quad_result = quadratize(poly_system, input_free=input_free, input_der_orders=input_der_orders,
                             conditions=conditions,
                             calc_upper_bound=calc_upper_bound,
                             generation_strategy=generation_strategy, scoring=scoring,
                             pruning_functions=pruning_functions,
                             new_vars_name=new_vars_name,
                             start_new_vars_with=start_new_vars_with + len(poly_system) - len(system))
    if quad_result:
        quad_result.polynomialization = poly_system
    return quad_result


def quadratize(poly_system: list[PolyElement] | list[(sp.Symbol, sp.Expr)] | EquationSystem,
               input_der_orders=None,
               input_free=False,
               conditions: Collection["SystemCondition"] = (),
               calc_upper_bound=True,
               generation_strategy=default_generation,
               scoring: Scoring = default_scoring,
               pruning_functions: Collection["Pruning"] | None = None,
               new_vars_name='w', start_new_vars_with=0) -> QuadratizationResult | None:
    """
    Quadratize a system of ODEs with the polynomial right-hand side.

    :param poly_system: system of polynomial equations in form [(X, p(X)), ...] where the left-hand side is derivatives of X.
    :param input_free: if True the function will not introduce derivatives of input functions.
    :param input_der_orders: a mapping of input variables to maximum order of their derivatives.
     For example, {T: 2} => T' and T'' are introduced.
    :param conditions: a list of predicates PolynomialSystem -> bool the quadratized systems must comply with.
     For example, use `partial(without_variables, [x1, x2])` to stop the algorithm from using x1 and x2 in new varaibles.
    :param calc_upper_bound: if True a non-optimal quadratization will be quickly found and used as upper bound.
     Disable if you already know an upper bound and add it in pruning_functions via `pruning_by_vars_number`.
    :param generation_strategy: a function for proposing new variables during quadratization evaluation.
    :param scoring: an ordering function for new variables.
    :param pruning_functions: a predicate indicating should the algorithm drop a current branch.
     Unlike `conditions` parameter, prunings apply to each intermediate result, not just to quadratizations.
     Check functions starting with `pruning_by` to see examples.
    :param new_vars_name: base name for new variables. For example, new_var_name='z' => z0, z1, z2, ...
    :param start_new_vars_with: initial index for new variables. Example: start_new_vars_with=3 => w_3, w_4, ...
    :return: a container of a quadratized system and new variables introduced or None if there is nothing found

    Example:
        >>> from qbee import *
        >>> x1, x2, u = functions("x1, x2, u")
        >>> quadratize(([x1, x1 + x1 * u), (x2, x1**2 * u)], input_free=True).print()
        Introduced variables:
        w0 = x1**2
         ‎
        x1' = u*x1 + x1
        x2' = u*w0
        w0' = 2*u*w0 + 2*w0
    """
    if pruning_functions is None:
        pruning_functions = default_pruning_rules

    if isinstance(poly_system, list) and isinstance(poly_system[0], tuple):
        poly_system = eq_list_to_eq_system(poly_system)

    if isinstance(poly_system, EquationSystem):
        if not poly_system.is_polynomial():
            raise Exception("Nonpolynomial system is passed to `quadratize` function.")

        if input_der_orders is None:
            if input_free:
                input_der_orders = {var: 0 for var in poly_system.variables.input if "'" not in str_qbee(var)}
            else:
                input_der_orders = {var: 1 for var in poly_system.variables.input if "'" not in str_qbee(var)}

        poly_equations, excl_inputs, all_inputs = poly_system.to_poly_equations(input_der_orders)
        without_excl_inputs = partial(without_variables, excl_vars=excl_inputs)
        pruning_by_decl_inputs = partial(pruning_by_declining_variables, excl_vars=excl_inputs)
        pruning_functions = [pruning_by_decl_inputs, *pruning_functions]
        conditions = [without_excl_inputs, *conditions]
    elif isinstance(poly_system, list) and isinstance(poly_system[0], PolyElement):
        poly_equations = poly_system
        excl_inputs, all_inputs = None, None
    else:
        raise TypeError("Incorrect type of the `system` parameter in quadratization.")

    result = quadratize_poly(poly_equations,
                             conditions=conditions,
                             calc_upper_bound=calc_upper_bound,
                             generation_strategy=partial(generation_strategy, excl_vars=excl_inputs),
                             scoring=scoring,
                             pruning_functions=[pruning_by_best_nvars, *pruning_functions],
                             new_vars_name=new_vars_name,
                             start_new_vars_with=start_new_vars_with)
    if result:
        result.exclude_variables(all_inputs)
    return result


def quadratize_poly(polynomials: list[PolyElement],
                    conditions: Collection["SystemCondition"] = (),
                    calc_upper_bound=True,
                    generation_strategy=default_generation,
                    scoring: Scoring = default_scoring,
                    pruning_functions: Collection["Pruning"] | None = None,
                    new_vars_name='w',
                    start_new_vars_with=0) -> QuadratizationResult | None:
    """
    Quadratize a system of ODEs with the polynomial right-hand side.

    :param polynomials: right-hand side of a system of polynomial ODEs, ordereded as in the variables in polynomials' Ring.
     Example: Ring("x, y, z", ...) => x' = p_0, y' = p_1, z' = p_2
    :param conditions: a list of predicates PolynomialSystem -> bool the quadratized systems must comply with.
     For example, use `partial(without_variables, [x1, x2])` to stop the algorithm from using x1 and x2 in new varaibles.
    :param calc_upper_bound: if True a non-optimal quadratization will be quickly found and used as upper bound.
     Disable if you already know an upper bound and add it in pruning_functions via `pruning_by_vars_number`.
    :param generation_strategy: a function for proposing new variables during quadratization evaluation.
    :param scoring: an ordering function for new variables.
    :param pruning_functions: a predicate indicating should the algorithm drop a current branch.
     Unlike `conditions` parameter, prunings apply to each intermediate result, not just to quadratizations.
     Check functions starting with `pruning_by` to see examples.
    :param new_vars_name: base name for new variables. For example, new_var_name='z' => z0, z1, z2, ...
    :param start_new_vars_with: initial index for new variables. Example: start_new_vars_with=3 => w_3, w_4, ...
    :return: a container of a quadratized system and new variables introduced or None if there is nothing found

    Example:
        >>> from sympy import ring, QQ
        >>> R, x, y = ring("x, y", QQ)
        >>> quad_res = quadratize_poly([x**2 * y, x * y**3])
        >>> quad_res.print()
        Introduced variables:
        w0 = x*y
        w1 = x*y**2

        x' = w0*x
        y' = w1*y
        w0' = w0**2 + w0*w1
        w1' = w0*w1 + 2*w1**2
    """
    if pruning_functions is None:
        pruning_functions = default_pruning_rules

    system = PolynomialSystem(polynomials)
    algo = BranchAndBound(system, conditions, generation_strategy, scoring,
                          (pruning_by_best_nvars,) + tuple(pruning_functions))
    if calc_upper_bound:
        algo.domination_upper_bound()
    algo_res = algo.quadratize()
    if algo_res.system is not None:
        quad_eqs, eq_vars, quad_vars = apply_quadratization(polynomials, algo_res.system.introduced_vars,
                                                            new_vars_name, start_new_vars_with)
        return QuadratizationResult(quad_eqs, eq_vars, quad_vars, algo_res)
    return None


# ------------------------------------------------------------------------------

class PolynomialSystem:
    def __init__(self, polynomials: list[PolyElement]):
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
            new_var_name + ("%d" % i) + " = " + monom2str(m, self.gen_symbols)
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
    def __init__(self, system: PolynomialSystem | None, introduced_vars: int, nodes_traversed: int):
        self.system = system
        self.num_introduced_vars = introduced_vars
        self.nodes_traversed = nodes_traversed

    def make_report(self, new_var_name="z_", start_new_vars_with=0) -> str:
        if self.system is None:
            return "No quadratization found under the given condition\n" + \
                   f"Nodes traversed: {self.nodes_traversed}"
        return f"Number of introduced variables: {self.num_introduced_vars}\n" + \
               f"Nodes traversed: {self.nodes_traversed}\n" + \
               "Introduced variables:\n" + self.system.to_str(new_var_name, start_new_vars_with)

    def __repr__(self):
        return self.make_report()


class QuadratizationResult:
    def __init__(self, equations, variables, quad_variables, quad_res: AlgorithmResult,
                 poly_res: EquationSystem | None = None):
        self.nodes_traversed = quad_res.nodes_traversed
        self.variables = variables
        self.equations = [sp.Eq(lhs, rhs.as_expr(), evaluate=False) for lhs, rhs in
                          zip(derivatives(variables), equations)]
        self.quadratization = quad_res.system
        self.polynomialization = poly_res
        self._excl_ders = []
        self._quad_variables = quad_variables

    @property
    def introduced_variables(self) -> list:
        poly_vars = self.polynomialization.introduced_variables if self.polynomialization else []
        quad_rhs = [monom2PolyElem(v, self.quadratization.gen_symbols) for v in self.quadratization.introduced_vars]
        quad_vars = [sp.Eq(dx, fx) for dx, fx in zip(self._quad_variables, quad_rhs)]
        return poly_vars + quad_vars

    def exclude_variables(self, variables):
        if variables and len(variables) > 0:
            self._excl_ders.extend(derivatives(variables))

    def print(self, str_function=str_qbee, with_introduced_variables=True):
        intr_vars_str = "Introduced variables:\n" + '\n'.join([str_function(eq) for eq in self.introduced_variables])
        equations_str = '\n'.join([
            str_function(eq) for eq in self.equations if eq.lhs not in self._excl_ders
        ])
        if with_introduced_variables:
            print(intr_vars_str + '\n\n' + equations_str)
        else:
            print(equations_str)

    def __repr__(self):
        intr_vars_str = '\n'.join([str_qbee(eq) for eq in self.introduced_variables])
        equations_str = '\n'.join([
            str_qbee(eq) for eq in self.equations if eq.lhs not in self._excl_ders
        ])
        return intr_vars_str + '\n\n' + equations_str

    @property
    def new_vars_count(self):
        return len(self.introduced_variables)


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

    @progress_bar(is_stop=False)
    def _dls(self, part_res: PolynomialSystem, to_depth: int, pred: Callable[[PolynomialSystem], bool], res: set):
        if part_res.new_vars_count() > to_depth:
            return

        if pred(part_res):
            res.add(part_res)
        else:
            for next_system in self.next_gen(part_res):
                self._dls(next_system, to_depth, pred, res)
        return

    @dump_results
    def get_optimal_quadratizations(self) -> Set[PolynomialSystem]:
        optimal_first = self.quadratize()
        return self.get_quadratizations(optimal_first.num_introduced_vars)

    @dump_results
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

    @logged(is_stop=False)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.generation_strategy.self.scoring)

    @progress_bar(is_stop=True)
    @logged(is_stop=True)
    def _final_iter(self):
        pass


# ------------------------------------------------------------------------------

QUAD_ALGORITHM_INTERRUPTED = False


def signal_handler(sig_num, frame):
    global POLY_ALGORITHM_INTERRUPTED
    print("The algorithm has been interrupted. Returning the current best.")
    QUAD_ALGORITHM_INTERRUPTED = True


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

    def quadratize(self) -> AlgorithmResult:
        nvars, opt_system, traversed = self._bnb_step(self._system, self.preliminary_upper_bound)
        self._final_iter()
        self._save_results(opt_system)
        return AlgorithmResult(opt_system, nvars, traversed)

    @progress_bar(is_stop=False)
    def _bnb_step(self, part_res: PolynomialSystem, best_nvars) \
            -> Tuple[int | float, PolynomialSystem | None, int]:
        self._nodes_traversed += 1
        # The order of these blocks is important: pruning rules assume that
        # the input partial result is not a quadratization
        if part_res.is_quadratized() and all(cond(part_res) for cond in self._sys_cond):
            return part_res.new_vars_count(), part_res, 1
        if any(map(lambda f: f(self, part_res, best_nvars), self._pruning_funs)) or QUAD_ALGORITHM_INTERRUPTED:
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

    @logged(is_stop=False)
    def next_gen(self, part_res: PolynomialSystem):
        return part_res.next_generation(self.generation_strategy, self.scoring)

    @progress_bar(is_stop=True)
    @logged(is_stop=True)
    def _final_iter(self):
        self._nodes_traversed = 0

    @dump_results
    def _save_results(self, opt_system):
        return [opt_system, ]


# ------------------------------------------------------------------------------

def pruning_by_nodes_processed(algo: Algorithm, _: PolynomialSystem, *args, nodes_processed: int):
    """
    Stops a search when the algorithm checks 'nodes_processed' nodes.

    Example:
        >>> from functools import partial
        >>> pruning = partial(pruning_by_nodes_processed, nodes_processed=100000)
    """
    if algo._nodes_traversed >= nodes_processed:
        return True
    return False


def pruning_by_nodes_without_quadratization_found(algo: Algorithm, _: PolynomialSystem, *args, nodes_processed: int):
    """
    Stops a search when the algorithm can not find a quadratization after checking `nodes_processed` nodes.

    Example:
        >>> from functools import partial
        >>> pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    """
    best_nvars = args[0]
    if best_nvars < math.inf:
        return False
    elif algo._nodes_traversed >= nodes_processed:
        return True
    return False


def pruning_by_elapsed_time(algo: Algorithm, system: PolynomialSystem, *args, start_t, max_t):
    """
    Stops a search after 'max_t' seconds.

    Example:
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

    Example:
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
    """Internal optimization pruning rule. For details check Section 5.1 from https://arxiv.org/abs/2103.08013"""
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


def pruning_by_declining_variables(a: Algorithm, part_res: PolynomialSystem, *args, excl_vars: list[Tuple]):
    """
    Prune out systems with `excl_vars` in quadratization

    Example:
        Tupples will correspond to variables ordering: z in Ring('x, y, z') => (0, 0, 1)

        >>> from functools import partial
        >>> pruning = partial(pruning_by_declining_variables, excl_vars=[(0, 0, 1)])
    """
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
    """Internal optimization pruning rule. For details check Section 5.2 from https://arxiv.org/abs/2103.08013"""
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
    """Internal optimization pruning rule that is used for a fast search of suboptimal quadratization."""
    if not part_res.introduced_vars:
        return False

    return not all([dominated(m, dominators) for m in part_res.introduced_vars])


default_pruning_rules = [
    pruning_by_quadratic_upper_bound,
    pruning_by_squarefree_graphs
]


# ------------------------------------------------------------------------------------------------

def without_variables(part_res: PolynomialSystem, excl_vars: list[Tuple]):
    """Deny quadratizations which have `excl_vars`."""
    excl_indices = [np.argmax(v) for v in excl_vars]
    if any([any([v[i] != 0 for i in excl_indices]) for v in part_res.introduced_vars]):
        return False
    return True


def with_higher_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) > system.original_degree, system.vars))


def with_le_degree_than_original(system: PolynomialSystem) -> bool:
    return any(map(lambda m: monomial_deg(m) <= system.original_degree, system.vars))
