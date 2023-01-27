from sympy import *
from qbee import *

if __name__ == '__main__':
    N = 6
    vs = functions(",".join([f"v_{i}" for i in range(1, N + 1)]))
    graph = {1: {2, 4}, 2: {3}, 3: {}, 4: {}, 5: {1, 2, 3, 6}, 6: {1, 2, 3}}
    eqs = []
    for i in range(1, N + 1):
        Av = sum([vs[j - 1] for j in graph[i]] + [vs[i - 1]])
        eqs.append([vs[i - 1], Av * vs[i - 1]**(-1) - Av])

    quadr_system = polynomialize_and_quadratize(eqs, generation_strategy=generation_semidiscretized)
    if quadr_system:
        print(quadr_system)
