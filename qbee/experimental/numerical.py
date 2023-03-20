from __future__ import annotations

import sympy as sp
import re
from examples import *
from functools import partial
from typing import Sequence
from ..quadratization import QuadratizationResult
from ..polynomialization import eq_list_to_eq_system
from ..printer import str_qbee
from ..util import tuple_to_monom
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
    params_sym = {sp.Symbol(str_qbee(k)): v for k, v in params.items()}

    if inputs is None:
        inputs = dict()
    inputs_sym = {sp.Symbol(str_qbee(k)): v for k, v in inputs.items()}

    system = eq_list_to_eq_system(res)
    states_sym = {sp.Symbol(str_qbee(k)): v for k, v in state_values.items()}

    def func(y, t, *args):
        vars_subs = dict(zip(states_sym.keys(), y))
        inputs_subs = {k: v.evalf(subs={INDEPENDENT_VARIABLE: t}) for k, v in inputs_sym.items()}
        return [eq.rhs.evalf(subs={**vars_subs, **params_sym, **inputs_subs}) for eq in system.equations]

    return partial(odeint, func, list(states_sym.values()))


def to_odeint_quadr(res: QuadratizationResult,
                    state_values: dict,
                    params: dict | None = None,
                    inputs: dict | None = None):
    if params is None:
        params = dict()
    param_subs = {sp.Symbol(str_qbee(k)): v for k, v in params.items()}

    if inputs is None:
        inputs = dict()
    inputs_sym = {sp.Symbol(str_qbee(k)): v for k, v in inputs.items()}
    inputs_ders_sym = generate_inputs_derivatives_subs(res, inputs_sym)

    states_sym = {sp.Symbol(str_qbee(k)): v for k, v in state_values.items()}
    poly_subs = {k: v.subs({**states_sym, **param_subs}) for k, v in res.polynomialization._substitution_equations.items()} \
        if res.polynomialization else dict()

    quad_vars_lhs = res._quad_variables
    quad_vars_rhs = [tuple_to_monom(m, res.quadratization.gen_symbols) for m in res.quadratization.introduced_vars]
    # Backward compatibility to Python 3.8 and less. Should be states_sym | poly_subs
    quad_subs = {k: v.subs({**states_sym, **poly_subs}) for k, v in zip(quad_vars_lhs, quad_vars_rhs)}

    # state_subs = states_sym | poly_subs | quad_subs
    state_subs = {**states_sym, **poly_subs, **quad_subs}

    def eval_input_subs_at(t):
        return {k: v.evalf(subs={INDEPENDENT_VARIABLE: t}) for k, v in {**inputs_sym, **inputs_ders_sym}.items()}
    def func(y, t, *args):
        vars_subs = dict(zip(state_subs.keys(), y))
        inputs_subs = eval_input_subs_at(t)
        return [eq.rhs.as_expr().evalf(subs={**vars_subs, **param_subs, **inputs_subs}) for eq in res.equations
                if eq.lhs not in res._excl_ders]

    inputs_at_0 = eval_input_subs_at(0)
    return partial(odeint, func, [v.evalf(subs={**inputs_at_0, INDEPENDENT_VARIABLE: 0.0})
                                  if isinstance(v, sp.Expr) else v
                                  for v in state_subs.values()])


def generate_inputs_derivatives_subs(res: QuadratizationResult, inputs_sym: dict[sp.Symbol, sp.Expr]) -> dict:
    lhs = [eq.lhs for eq in res.equations]
    inputs_ders = [apply_derivative(subs, get_input_derivatives(lhs, inp)) for inp, subs in inputs_sym.items()]
    return reduce(lambda a, b: {**a, **b}, inputs_ders, dict())

def apply_derivative(inp_subs: sp.Expr, ders: list[sp.Symbol]) -> dict[sp.Symbol, sp.Expr]:
    return {du: inp_subs.diff(INDEPENDENT_VARIABLE, str_qbee(du).count('\'')) for du in ders}


def get_input_derivatives(symbols: list[sp.Symbol], inp: sp.Symbol):
    return [s for s in symbols if re.match(rf"^{inp}\'+$", str_qbee(s))]
