import sympy as sp
from qbee import *
from benchmark.main import benchmark_system

x = sp.symbols('x')
x_dot = derivatives(x)[0]

system = EquationSystem([
    sp.Eq(x_dot, x / (1 + sp.exp(x)), evaluate=False)
])
system = PolynomialSystem.from_EquationSystem(polynomialize(system))


def xSigmoid_benchmark():
    restricted_heuristics = ('none', 'FVC', 'AED', 'AEQD', 'SMD')
    benchmark_system(system=system, system_name="xSigmoid", cycles=20, search_algorithms=('BFS', 'ID-DLS'),
                     heuristics=restricted_heuristics, limit_depth=2)


if __name__ == '__main__':
    xSigmoid_benchmark()
