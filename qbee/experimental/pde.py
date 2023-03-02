from __future__ import annotations

import sympy as sp
from sympy.core.function import AppliedUndef
from copy import deepcopy
from queue import Queue
from typing import Optional, Collection
from functools import reduce
from qbee.printer import str_qbee
from qbee.quadratization import SystemCondition, QuadratizationResult, polynomialize_and_quadratize, Pruning
from qbee.selection import default_scoring, default_generation, Scoring


def polynomialize_and_quadratize_pde(start_system: list[(sp.Symbol, sp.Expr)],
                                     input_der_orders: dict | None = None,
                                     conditions: Collection["SystemCondition"] = (),
                                     polynomialization_upper_bound=10,
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
        # if pb_enable:
        #     print("Current spatial time derivatives equations:")
        #     print("...")
        #     for eq in system[len(start_system):]:
        #         print_qbee(eq)
        #     print()

        quad_res = polynomialize_and_quadratize(system, input_orders_with_pde, conditions,
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


def select_inputs(system):
    lhs, rhs = zip(*system)
    funcs = set(reduce(lambda l, r: l | r, [eq.atoms(AppliedUndef) for eq in rhs]))
    lhs_args = set(sp.flatten([eq.args for eq in lhs if not isinstance(eq, sp.Derivative)]))
    return set(filter(lambda f: (f not in lhs) and (f not in lhs_args), funcs))


def select_pde_inputs(system):
    lhs, rhs = zip(*system)
    return set(filter(lambda v: v not in lhs, reduce(lambda l, r: l | r, [eq.atoms(sp.Derivative) for eq in rhs])))


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
