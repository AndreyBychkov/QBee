# Computes the largest number of edges in C4*-free graphs, 
# see Table 1 in the paper (https://arxiv.org/pdf/2103.08013.pdf)

import math
from itertools import combinations


def has_C4(G):
    common_heighbours = set()
    for neighbours in G:
        for i in range(len(neighbours)):
            for j in range(i + 1, len(neighbours)):
                u = neighbours[i]
                v = neighbours[j]
                pair = (min(u, v), max(u, v))
                if pair in common_heighbours:
                    return True
                common_heighbours.add(pair)
    return False


def graph_from_edges(n, edges):
    G = [[] for _ in range(n)]
    for i, j in edges:
        G[i].append(j)
        if i != j:
            G[j].append(i)
    return G


def max_num_edges(n):
    result = []

    # without loss of generality, we assume that the first two vertices are connected
    edges = [(0, 1)]
    edges_to_choose = []
    for i in range(n):
        for j in range(i + 1, n):
            if i != 0 or j != 1:
                edges_to_choose.append((i, j))

    for max_num_loops in range(n):
        best_so_far, witness = 0, []
        if len(result) > 0:
            best_so_far, witness = result[-1]
        # without loss of generality, we assume that the loops are on the nodes
        # 1, ..., max_num_loops
        if max_num_loops != 0:
            edges.append((max_num_loops, max_num_loops))
        taken = len(edges)

        for to_add in range(max(0, best_so_far - taken + 1), len(edges_to_choose) + 1):
            found_something = False
            for new_edges in combinations(edges_to_choose, to_add):
                G = graph_from_edges(n, edges + list(new_edges))
                if not has_C4(G):
                    best_so_far = to_add + taken
                    witness = edges + list(new_edges)
                    found_something = True
                    break
            if not found_something:
                break
        result.append((best_so_far, witness))

    return result


####################

for i in range(2, 8):
    print(i)
    print("\n".join(map(str, max_num_edges(i))))
