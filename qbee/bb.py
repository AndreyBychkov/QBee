import copy
import functools as ft
import math
import time

from sympy import *

def get_decompositions(monomial):
    if len(monomial) == 0:
        return {(tuple(), tuple())}
    result = set()
    prev_result = get_decompositions(tuple(monomial[:-1]))
    for r in prev_result:
        for i in range(monomial[-1] + 1):
            a, b = tuple(list(r[0]) + [i]), tuple(list(r[1]) + [monomial[-1] - i])
            result.add((min(a, b), max(a, b)))
    return result
        
class PartialResult:
    def __init__(self, polynomials):
        """
        polynomials - right-hand sides of the ODE system listed in the same order as 
                      the variables in the polynomial ring
        """
        self.dim = len(polynomials[0].ring.gens)

        # put not monomials but differents in the exponents between rhs and lhs
        self.rhs = dict()
        for i, p in enumerate(polynomials):
            self.rhs[i] = set()
            for m in p.to_dict().keys():
                mlist = list(m)
                mlist[i] -= 1
                self.rhs[i].add(tuple(mlist))

        self.vars, self.squares, self.nonsquares = set(), set(), set()
        self.add_var(tuple([0] * self.dim))
        for i in range(self.dim):
            self.add_var(tuple([1 if i == j else 0 for j in range(self.dim)]))

    def add_var(self, v):
        for i in range(self.dim):
            if v[i] > 0:
                for m in self.rhs[i]:
                    self.nonsquares.add(tuple([v[j] + m[j] for j in range(self.dim)]))

        self.vars.add(v)
        for u in self.vars:
            self.squares.add(tuple([u[i] + v[i] for i in range(self.dim)]))
        self.nonsquares = set(filter(lambda s: s not in self.squares, self.nonsquares))

    def is_quadratized(self):
        return not self.nonsquares

    def get_system_score(self):
        total_nonsquare = sum([sum(m) for m in self.nonsquares])
        return total_nonsquare + self.dim * len(self.vars)

    def get_smallest_nonsquare(self):
        return min([(sum(m), m) for m in self.nonsquares])[1]

    def next_generation(self):
        new_gen = []
        for d in get_decompositions(self.get_smallest_nonsquare()):
            c = copy.deepcopy(self)
            for v in d:
                c.add_var(v)
            new_gen.append(c)

        return sorted(new_gen, key=lambda s: s.get_system_score())

    def new_vars_count(self):
        return len(self.vars) - self.dim  - 1

def find_optimal_quadratization(part_res, best=math.inf):
    if part_res.is_quadratized():
        return (part_res.new_vars_count(), part_res, 1)
    if part_res.new_vars_count() >= best - 1:
        return (math.inf, None, 1)

    traversed_total = 1
    min_nvars, best_system = best, None
    for next_system in part_res.next_generation():
        nvars, opt_system, traversed = find_optimal_quadratization(next_system, min_nvars)
        traversed_total += traversed
        if nvars < min_nvars:
            min_nvars = nvars
            best_system = opt_system
    return (min_nvars, best_system, traversed_total)

if __name__ == "__main__":
    R, x, y = ring(["x", "y"], QQ)

    system = PartialResult([x**2 + y, 3 * x * y**3 - y])
    start = time.time()
    nvars, sys, nodes = find_optimal_quadratization(system)
    end = time.time()
    print(f"Runtime: {end - start}")
    print(f"Order of quadratization {nvars}")
    print(f"All the variables {sys.vars}")
    print(f"Number of nodes traversed by the algo {nodes}")
    print("===========")

    system = PartialResult([(x + 1)**8, y])
    start = time.time()
    nvars, sys, nodes = find_optimal_quadratization(system)
    end = time.time()
    print(nvars)
    print(sys.vars)
    print(nodes)
    print(f"Runtime: {end - start}")
    print(f"Order of quadratization {nvars}")
    print(f"All the variables {sys.vars}")
    print(f"Number of nodes traversed by the algo {nodes}")
    print("===========")

    system = PartialResult([(x + 1)**15, y])
    start = time.time()
    nvars, sys, nodes = find_optimal_quadratization(system)
    end = time.time()
    print(nvars)
    print(sys.vars)
    print(nodes)
    print(f"Runtime: {end - start}")
    print(f"Order of quadratization {nvars}")
    print(f"All the variables {sys.vars}")
    print(f"Number of nodes traversed by the algo {nodes}")
    print("===========")

    # xSigmoid
    R, x, y, z = ring(["x", "y", "z"], QQ)
    system = PartialResult([x * z, x * y * z, x * y * z**3])
    start = time.time()
    nvars, sys, nodes = find_optimal_quadratization(system)
    end = time.time()
    print(nvars)
    print(sys.vars)
    print(nodes)
    print(f"Runtime: {end - start}")
    print(f"Order of quadratization {nvars}")
    print(f"All the variables {sys.vars}")
    print(f"Number of nodes traversed by the algo {nodes}")
    print("===========")

    # Rabinovich-Fabrikant
    R, x, y, z = ring(["x", "y", "z"], QQ)
    system = PartialResult([y * (z - 1 - x**2) + x, x * (3 * z + 1 - x**2) + y, -2 * z * (2 + x * y)])
    start = time.time()
    nvars, sys, nodes = find_optimal_quadratization(system)
    end = time.time()
    print(nvars)
    print(sys.vars)
    print(nodes)
    print(f"Runtime: {end - start}")
    print(f"Order of quadratization {nvars}")
    print(f"All the variables {sys.vars}")
    print(f"Number of nodes traversed by the algo {nodes}")
    print("===========")   
 
