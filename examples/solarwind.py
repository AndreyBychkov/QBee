from qbee import *

if __name__ == '__main__':
    N = 4
    vs = functions(",".join([f"v{i}" for i in range(N)]))
    us = functions(",".join([f"u{i}" for i in range(N)]))
    graph = {0: {0, 1, 3}, 1: {1, 2}, 2: {2}, 3: {3}}
    system = []
    for i in range(N):
        system.append((vs[i], sum([us[j] for j in graph[i]]) + sum([vs[j] for j in graph[i]])))
        system.append((us[i], system[-1][1] / vs[i]))

    quadr_system = polynomialize_and_quadratize(system)
    if quadr_system:
        print(quadr_system)
