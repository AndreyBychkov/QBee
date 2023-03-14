from sympy import *
from qbee import *

from functools import partial

if __name__ == '__main__':

    n = 6
    R, *allvars = ring(",".join([f"ps_{i}, th_{i}, w_{i}" for i in range(1, n + 1)]) + ", u, up", QQ)
    
    def psi(i):
        return allvars[(i - 1) * 3]

    def theta(i):
        return allvars[(i - 1) * 3 + 1]

    def w(i):
        return allvars[(i - 1) * 3 + 2]

    u, up = allvars[-2:]

    graph = {1: {2, 3}, 2: set(), 3: set(), 4: {5, 6}, 5: {4, 6}, 6: {4, 5}}
    eqs = []
    for i in range(1, n + 1):
        eqs.extend(
            [psi(i) + u + sum([psi(j) for j in graph[i]]) + psi(i) * w(i),
            theta(i) + u + sum([theta(j) for j in graph[i]]) + psi(i) * w(i)]
        )
        eqs.append(eqs[-1] * theta(i)**(-2) * w(i))

    eqs.extend([up, up])
    up_as_tuple = tuple([0] * (len(allvars) - 1) + [1])
    quadr_system = quadratize(eqs, generation_strategy=partial(generation_semidiscretized, excl_vars=[up_as_tuple]))
    if quadr_system:
        quadr_system.print()
