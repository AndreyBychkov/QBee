import sympy as sp
from qbee import *
from benchmark.main import benchmark_system

x, y, z, w = sp.symbols('x, y, z, w')
x_dot, y_dot, z_dot, w_dot = derivatives('x, y, z, w')

system = EquationSystem([
    sp.Eq(x_dot, y - x + y * z * w),
    sp.Eq(y_dot, x + z - x * z * w),
    sp.Eq(z_dot, -z + x * y * w),
    sp.Eq(w_dot, -w + x * y * z)
])

if __name__ == '__main__':
    benchmark_system(system=system, system_name="4DQiDuChenChenYuan", cycles=10, search_algorithms=('BFS', 'ID-DLS'),
                     heuristics=('none', 'FF', 'FVC', 'AED', 'AEQD', 'SMD'), initial_max_depth=4, limit_depth=4)