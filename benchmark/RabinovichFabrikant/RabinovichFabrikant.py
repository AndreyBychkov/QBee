import sympy as sp
from qbee import *
from benchmark.main import benchmark_system

x, y, z = sp.symbols('x, y, z')
x_dot, y_dot, z_dot = derivatives('x, y, z')

system = EquationSystem([
    sp.Eq(x_dot, y * (z - 1 + x ** 2) + x),
    sp.Eq(y_dot, x * (3 * z + 1 - x ** 2) + y),
    sp.Eq(z_dot, -2 * z * (1 + x * y))
])

if __name__ == '__main__':
    benchmark_system(system=system, system_name="RabinovichFabrikant", cycles=20, search_algorithms=('BFS', 'ID-DLS'),
                     heuristics=('none', 'FF', 'FVC', 'AED', 'AEQD', 'SMD'), initial_max_depth=3, limit_depth=3)