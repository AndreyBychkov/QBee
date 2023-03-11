from __future__ import annotations

import sympy as sp
from .util import *
from .quadratization import *

# The graph from Proposition 4.7
DA_GRAPH = {0: {0, 1, 3}, 1: {1, 2}, 2: {2}, 3: {3}}


def quadratize_dimension_agnostic(system: list[(sp.Symbol, sp.Expr)], 
                                  non_duplicated_vars: list[sp.Symbol] | None = None,
                                  input_der_orders=None,
                                  input_free=False,
                                  print_intermediate=False):
    """
    Find a dimension-agnostic quadratization for a linearly coupled family of ODE system

    :param system: the symbolic representation of the input family, where Dx is the coupling term for the variables x
    :param non_duplicated_vars: list of variables which are the same in all the blocks of the family (must be inputs)
    :param input_der_orders: maximal oders for the input variables (if defined, overwrites input_free)
    :param input_free: if an input-free quadratization is sought
    :param print_intermediate: if the quadratization of the auxiliary ODE system from Proposition 4.7 should be printed
    :return: dimension-agnistic quadratization of the input family
    """

    if non_duplicated_vars is None:
        non_duplicated_vars = []
    odesubs = [{} for _ in range(len(DA_GRAPH))]
    all_funcs = reduce(lambda l, r: l | r, map(lambda e: e.atoms(sp.Function), [e[0] for e in system] + [e[1] for e in system]))

    # duplicating the variables
    for f in all_funcs:
        if str_qbee(f)[0] != 'D' and f not in non_duplicated_vars:
            for i in range(len(DA_GRAPH)):
                odesubs[i][f] = functions(str_qbee(f) + "_" + str(i), laurent=f.is_laurent)

    # defining the coupling
    for Df in all_funcs:
        if str_qbee(Df)[0] == 'D':
             f = next(F for F in all_funcs if str_qbee(F) == str_qbee(Df)[1:])
             for i in range(len(DA_GRAPH)):
                 odesubs[i][Df] = sum([odesubs[j][f] for j in DA_GRAPH[i]])

    # constructing the system from Proposition 4.7
    ode_system = []
    for lhs, rhs in system:
        if lhs in non_duplicated_vars:
            error(f"Non-duplicated variables must not appear on the left-hand side, not true for {lhs}")
        for i in range(len(DA_GRAPH)):
            ode_system.append((lhs.subs(odesubs[i]), rhs.subs(odesubs[i])))

    # transforming the input info
    ode_input_der_orders = input_der_orders
    if input_der_orders is not None:
        ode_input_der_orders = {}
        for u, order in input_der_orders.items():
            if u in non_duplicated_vars:
                ode_input_der_orders[u] = order
            else:
                for i in range(len(DA_GRAPH)):
                    ode_input_der_orders[odesubs[i][u]] = order

    # computing ODE quadratization
    quadr_system = quadratize(ode_system, 
                              generation_strategy=generation_semidiscretized,
                              input_free=input_free,
                              input_der_orders=ode_input_der_orders)
    if print_intermediate:
        print("Quadratization of the auxiliary ODE system:")
        quadr_system.print()

    # extracting dimension-agnostic quadratization
    da_new_vars = []
    uncoupled = {sp.symbols(str_qbee(v)) : sp.symbols(str_qbee(k)) for k, v in odesubs[0].items()}
    coupled = {sp.symbols(str_qbee(v)) : sp.symbols(str_qbee(k) + "_tilda") for k, v in odesubs[1].items()}
    non_duplicated = {sp.symbols(str_qbee(v)) for v in non_duplicated_vars}
    for eq in quadr_system.introduced_variables:
        new_var = eq.rhs
        if new_var.free_symbols.issubset(uncoupled.keys() | non_duplicated):
            da_new_vars.append(new_var.subs(uncoupled))
        elif new_var.free_symbols.issubset(uncoupled.keys() | non_duplicated | coupled.keys()) and not new_var.free_symbols.issubset(non_duplicated | coupled.keys()):
            da_new_vars.append(new_var.subs(uncoupled | coupled))

    return da_new_vars
