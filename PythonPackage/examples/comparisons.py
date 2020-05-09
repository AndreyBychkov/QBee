import sympy as sp
from qbee import *

sp.init_printing()

x = sp.symbols('x')
dot_x = make_derivative_symbol(x)

system = EquationSystem([
    sp.Eq(dot_x, x / (1 + sp.exp(x)))
])

system = polynomialize(system)


def ql_bfs_stat(_) -> int:
    ql_system = quadratic_linearize(system, mode="optimal", method_optimal="bfs")
    return ql_system.statistics.steps


def ql_iddfs_1_stat(_) -> int:
    ql_system = quadratic_linearize(system, mode="optimal", method_optimal="iddfs",
                                    initial_max_depth=1, heuristics='replacement-value')
    return ql_system.statistics.steps


def ql_iddfs_2_stat(_) -> int:
    ql_system = quadratic_linearize(system, mode="optimal", method_optimal="iddfs",
                                    initial_max_depth=2, heuristics='replacement-value')
    return ql_system.statistics.steps
