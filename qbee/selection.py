from sympy.polys.monomials import monomial_deg, monomial_divides
from typing import Tuple, Callable, Iterable, Any

from .util import *

GenerationStrategy = Callable[['PolynomialSystem'], Iterable[Any]]
Scoring = Callable[['PolynomialSystem'], int]

# Standard generation strategy with different scoring functions

def default_generation(system):
    if len(system.nonsquares) == 0:
        return list()
    return get_decompositions(system.get_smallest_nonsquare())

#####

"""
Given two mappings from the variables to indices and a list of numbers
creates a new list in which the values at the indices corresponding to 
variables in `v_to_ind_from` are placed at positions as prescribed by `v_to_ind_to`
"""
def _map_indices(v_to_ind_from, v_to_ind_to, arr, no_class):
    result = [arr[i] if i in no_class else 0 for i in range(len(arr))]
    v_to_val = {v: arr[ind] for v, ind in v_to_ind_from.items()}
    for v, ind in v_to_ind_to.items():
        result[ind] = v_to_val[v]
    return result


def generation_semidiscretized(system):
    # recognizing the semidiscretized structure if not cashed already
    if not hasattr(system, "equivalence_classes"):
        system.equivalence_classes = dict()
        system.ind_to_class = dict()
        system.ind_to_var = dict()
        system.no_class = set()
        for i, s in enumerate(system.gen_symbols):
            if '_' not in str(s):
                system.no_class.add(i)
                continue
            class_ind = int(str(s).split("_")[-1])
            var_name = "_".join(str(s).split("_")[:-1])
            if class_ind not in system.equivalence_classes:
                system.equivalence_classes[class_ind] = dict()
            system.equivalence_classes[class_ind][var_name] = i
            system.ind_to_class[i] = class_ind
            system.ind_to_var[i] = var_name

        # building the graph
        system.graph = dict()
        for i, eq in system.rhs.items():
            if i in system.no_class:
                continue
            c = system.ind_to_class[i]
            if c not in system.graph:
                system.graph[c] = set()
            for m in eq:
                classes = set([system.ind_to_class[j] for j in system.ind_to_class.keys() if m[j] != 0])
                for cl in classes:
                    if cl != c:
                        system.graph[c].add(cl)

        print(system.equivalence_classes)
        print(system.ind_to_class)
        print(system.ind_to_var)
        print(system.graph)
    
    # producing sets of new variables
    decomposition = get_decompositions(system.get_smallest_nonsquare())
    result = set()
    for d in decomposition:
        extended_vars = set(d)
        for v in d:
            classes = set([system.ind_to_class[j] for j in system.ind_to_class.keys() if v[j] != 0])
            if len(classes) == 1:
                c = list(classes)[0]
                for cl, var_to_ind in system.equivalence_classes.items():
                    extended_vars.add(tuple(_map_indices(
                        system.equivalence_classes[c], var_to_ind, v, system.no_class
                    )))
            elif len(classes) == 2:
                c1, c2 = list(classes)
                for u1, adj in system.graph.items():
                    for u2 in adj:
                        for p in [(c1, c2), (c2, c1)]:
                            if p[1] in system.graph[p[0]]:
                                m1 = _map_indices(
                                    system.equivalence_classes[p[0]], system.equivalence_classes[u1], v, system.no_class
                                )
                                m2 = _map_indices(
                                    system.equivalence_classes[p[1]], system.equivalence_classes[u2], v, system.no_class
                                )
                                extended_vars.add(tuple([a + b for a, b in zip(m1, m2)]))
            elif len(classes) == 0:
                extended_vars.add(v)
            else:
                extended_vars = set()
                break
        if len(extended_vars) > 0:
             result.add(tuple(sorted(list(extended_vars))))
    return result


# Different scoring functions

def empty_score(system) -> int:
    return 1

def default_scoring(system) -> int:
    total_nonsquare = sum([sum(map(abs, m)) for m in system.nonsquares])
    return total_nonsquare + system.dim * len(system.vars)


def aeqd_scoring(system) -> int:
    eq_degs = list(map(lambda mlist: max(map(monomial_deg, mlist)), system.rhs.values()))
    aeqds = list(map(_compute_aeqd, system.nonsquares, [eq_degs, ] * len(system.nonsquares)))
    return sum(aeqds)


def smd_scoring(system) -> int:
    mlist = system.nonsquares
    return sum(map(lambda s: _compute_smd(s, mlist), mlist))


def _compute_scoring(sub: Tuple[int], eq_degs):
    mon_degs = map(lambda deg: deg + monomial_deg(sub) - 1, eq_degs)
    quad_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, mon_degs))
    return sum(quad_discrepancies)


def _compute_scoring(sub, mlist: list):
    return (monomial_deg(sub) - 1) * len(list(filter(lambda m: monomial_divides(sub, m), mlist)))
