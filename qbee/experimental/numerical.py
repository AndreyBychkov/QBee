from __future__ import annotations

import sympy as sp
from examples import *
from functools import partial
from typing import Sequence
from qbee.quadratization import QuadratizationResult
from qbee.polynomialization import eq_list_to_eq_system
from scipy.integrate import odeint


def to_odeint(system: Sequence[(sp.Symbol, sp.Expr)] | QuadratizationResult,
              state_values: dict,
              params: dict | None = None,
              inputs: dict | None = None):
    if isinstance(system, QuadratizationResult):
        return to_odeint_quadr(system, state_values, params, inputs)
    elif isinstance(system, Sequence) and isinstance(next(iter(system)), tuple):
        return to_odeint_no_quadr(system, state_values, params, inputs)
    else:
        raise TypeError("Incorrect type of the `system` parameter: "
                        "expected QuadratizationResult or Sequence[(sympy.Symbol, sympy.Expr)]")


def to_odeint_no_quadr(res: Sequence[(sp.Symbol, sp.Expr)],
                       state_values: dict,
                       params: dict | None = None,
                       inputs: dict | None = None):
    if params is None:
        params = dict()
    params_sym = {sp.Symbol(str(k)): v for k, v in params.items()}

    if inputs is None:
        inputs = dict()
    inputs_sym = {sp.Symbol(str(k)): v for k, v in inputs.items()}

    system = eq_list_to_eq_system(res)
    states_sym = {sp.Symbol(str(k)): v for k, v in state_values.items()}

    def func(y, t, *args):
        vars_subs = dict(zip(states_sym.keys(), y))
        inputs_subs = {k: v.evalf(subs={INDEPENDENT_VARIABLE: t}) for k, v in inputs_sym.items()}
        return [eq.rhs.evalf(subs=vars_subs | params_sym | inputs_subs) for eq in system.equations]

    return partial(odeint, func, list(states_sym.values()))


def to_odeint_quadr(res: QuadratizationResult,
                    state_values: dict,
                    params: dict | None = None,
                    inputs: dict | None = None):
    if params is None:
        params = dict()
    param_subs = {sp.Symbol(str(k)): v for k, v in params.items()}

    if inputs is None:
        inputs = dict()
    inputs_sym = {sp.Symbol(str(k)): v for k, v in inputs.items()}

    states_sym = {sp.Symbol(str(k)): v for k, v in state_values.items()}
    poly_subs = {k: v.subs(states_sym) for k, v in res.polynomialization._substitution_equations.items()}

    base_name = res.polynomialization.variables.base_var_name
    quad_start_index = res.polynomialization.variables.start_new_vars_with + \
                       len(res.polynomialization.variables.generated)
    quad_vars_lhs = sp.symbols(
        f"{base_name}{quad_start_index}:{quad_start_index + res.quadratization.new_vars_count()}")
    quad_vars_rhs = [tuple_to_monom(m, res.quadratization.gen_symbols) for m in res.quadratization.introduced_vars]
    quad_subs = {k: v.subs(states_sym) for k, v in zip(quad_vars_lhs, quad_vars_rhs)}

    state_subs = states_sym | poly_subs | quad_subs

    def func(y, t, *args):
        vars_subs = dict(zip(state_subs.keys(), y))
        inputs_subs = {k: v.evalf(subs={INDEPENDENT_VARIABLE: t}) for k, v in inputs_sym.items()}
        return [expr.as_expr().evalf(subs=vars_subs | param_subs | inputs_subs) for dx, expr in zip(res.lhs, res.rhs)
                if dx not in res.excl_ders]

    return partial(odeint, func, list(state_subs.values()))
