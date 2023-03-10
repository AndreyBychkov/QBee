from __future__ import annotations

import sympy as sp
from .util import *
from .quadratization import *

# The graph from Proposition 4.7
DA_GRAPH = {0: {0, 1, 3}, 1: {1, 2}, 2: {2}, 3: {3}}


def quadratize_dimension_agnostic(system: list[(sp.Symbol, sp.Expr)], non_duplicated_vars: list[sp.Symbol] | None = None):
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

    ode_system = []
    for lhs, rhs in system:
        if lhs in non_duplicated_vars:
            error(f"Non-duplicated variables must not appear on the left-hand side, not true for {lhs}")
        for i in range(len(DA_GRAPH)):
            ode_system.append((lhs.subs(odesubs[i]), rhs.subs(odesubs[i])))
   
    print(ode_system)
    quadr_system = quadratize(ode_system, generation_strategy=generation_semidiscretized)
    quadr_system.print()

