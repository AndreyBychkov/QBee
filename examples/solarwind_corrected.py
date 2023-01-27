from sympy import *
from qbee import *
from functools import partial

if __name__ == '__main__':
    N = 6
    R, *vs = ring(",".join([f"v_{i}, u_{i}" for i in range(1, N + 1)]), QQ)
    graph = {1: {2, 4}, 2: {3}, 3: {}, 4: {}, 5: {1, 2, 3, 6}, 6: {1, 2, 3}}
    def v(i):
        return vs[2 * i - 2]
    def u(i):
        return vs[2 * i - 1]

    eqs = []
    for i in range(1, N + 1):
        Dv = sum([v(j) for j in graph[i]] + [v(i)])
        Du = sum([u(j) for j in graph[i]] + [u(i)])
        eqs.append(Du - Dv)
        eqs.append(Du / v(i) - Dv / v(i))

    quadr_system = quadratize(eqs, generation_strategy=generation_semidiscretized)
    if quadr_system:
        print(quadr_system)
