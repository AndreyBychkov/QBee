from __future__ import annotations

import sympy as sp
from .util import *

# The graph from Proposition 4.7
DA_GRAPH = {0: {0, 1, 3}, 1: {1, 2}, 2: {2}, 3: {3}}


def quadratize_dimension_agnostic(system: list[(sp.Symbol, sp.Expr)]):
    main_vars = [eq[0] for eq in system]
    ode_system = []

    def ode_var(i, j):
        return sp.symbols(f"{main_vars[i]}_{j}")

    for i in range(len(DA_GRAPH)):
        sub = {v: ode_var(j, i) for j, v in enumerate(main_vars)}
        coupling = {sp.symbols(f"D{v}"): sum([ode_var(j, s) for s in DA_GRAPH[i]]) for j, v in enumerate(main_vars)}
        sub.update(coupling)
        for lhs, rhs in system:
            ode_system.append((lhs.subs(sub), rhs.subs(sub)))
    print(ode_system)

