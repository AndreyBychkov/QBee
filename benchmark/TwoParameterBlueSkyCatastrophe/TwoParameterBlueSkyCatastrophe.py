import sympy as sp
from qbee import *
from benchmark.main import benchmark_system

x, y, z = sp.symbols('x, y, z')
x_dot, y_dot, z_dot = derivatives('x, y, z')

system = EquationSystem([
    sp.Eq(x_dot, (2 - 10 * (x ** 2 + y ** 2)) * x + y ** 2 + 2 * y + z ** 2),
    sp.Eq(y_dot, -z ** 3 - (1 + y) * (y ** 2 + 2 * y + z ** 2) - 4 * x + y),
    sp.Eq(z_dot, (1 + y) * z ** 2 + x ** 2)
])

if __name__ == '__main__':
    all_heuristics = ('none', 'FF', 'FVC', 'AED', 'AEQD', 'SMD')
    restricted_heuristics = ('none', 'AED', 'AEQD', 'SMD')
    benchmark_system(system=system, system_name="TwoParameterBlueSkyCatastrophe", cycles=10, search_algorithms=('BFS', 'ID-DLS', 'MMDR'),
                     heuristics=restricted_heuristics, initial_max_depth=4, limit_depth=4)
