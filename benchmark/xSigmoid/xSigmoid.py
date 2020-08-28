import sympy as sp
from qbee import *
from benchmark.main import benchmark_system

x = sp.symbols('x')
x_dot = derivatives(x)

system = EquationSystem([
    sp.Eq(x_dot, x / (1 + sp.exp(x)))
])
system = polynomialize(system)


def xSigmoid_benchmark():
    restricted_heuristics = ('none', 'FVC', 'AED', 'AEQD', 'SMD')
    benchmark_system(system=system, system_name="xSigmoid", cycles=20, search_algorithms=('BFS', 'ID-DLS'),
                     heuristics=restricted_heuristics, initial_max_depth=2, limit_depth=2)


if __name__ == '__main__':
    xSigmoid_benchmark()
