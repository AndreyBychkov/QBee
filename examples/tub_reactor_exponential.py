from sympy import *
from qbee import *

if __name__ == '__main__':
    
    R, *allvars = ring(",".join([f"ps_{i}, th_{i}, w_{i}" for i in range(1, 4)]), QQ)
    def psi(i):
        return allvars[(i - 1) * 3]

    def theta(i):
        return allvars[(i - 1) * 3 + 1]

    def w(i):
        return allvars[(i - 1) * 3 + 2]


    graph = {1: {2, 3}, 2: set(), 3: set()}
    eqs = []
    for i in range(1, 4):
        eqs.extend(
            [psi(i) + sum([psi(j) for j in graph[i]]) + psi(i) * w(i),
            theta(i) + sum([theta(j) for j in graph[i]]) + psi(i) * w(i)]
        )
        eqs.append(eqs[-1] * theta(i)**(-2) * w(i))

    quadr_system = quadratize(eqs, generation_strategy=generation_semidiscretized)
    if quadr_system:
        print(quadr_system)
