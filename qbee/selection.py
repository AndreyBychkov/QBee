import numpy as np
from sympy.polys.monomials import monomial_deg, monomial_divides
from typing import Tuple, Callable, Iterable, Any

from .util import *

Scoring = Callable[['PolynomialSystem'], int]


def raw_generation(system, excl_vars=None):
    """Function for proposing some new variables (to be used in both generators"""
    if excl_vars is None:
        excl_vars = list()
    excl_indices = [np.argmax(v) for v in excl_vars]
    all_decomp = get_decompositions(system.get_smallest_nonsquare())
    if len(excl_indices) == 0:
        return all_decomp
    result = []
    for d in all_decomp:
        to_add = True
        for m in d:
            if sum(m) != 1 or sum(map(abs, m)) != 1:
                if len([i for i in excl_indices if m[i] != 0]) > 0:
                    to_add = False
                    break
        if to_add:
            result.append(d)
    return result


def default_generation(system, excl_vars=None):
    """Standard generation strategy with different scoring functions"""
    if excl_vars is None:
        excl_vars = list()
    if len(system.nonsquares) == 0:
        return list()
    return raw_generation(system, excl_vars)


def _map_indices(v_to_ind_from, v_to_ind_to, arr, no_class):
    """
    Given two mappings from the variables to indices and a list of numbers
    creates a new list in which the values at the indices corresponding to
    variables in `v_to_ind_from` are placed at positions as prescribed by `v_to_ind_to`
    """
    result = [arr[i] if i in no_class else 0 for i in range(len(arr))]
    v_to_val = {v: arr[ind] for v, ind in v_to_ind_from.items()}
    for v, ind in v_to_ind_to.items():
        result[ind] = v_to_val[v]
    return result


def generation_semidiscretized(system, excl_vars=None):
    if excl_vars is None:
        excl_vars = list()
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
            # stripping \' to handle derivatives of inputs
            var_split = str(s).split("_")
            class_ind = int(var_split[-1].strip('\''))
            var_name = "_".join(var_split[:-1]) + '\'' * var_split[-1].count('\'')
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

    # producing sets of new variables
    decomposition = raw_generation(system, excl_vars)
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


def _compute_aeqd(sub: Tuple[int], eq_degs):
    mon_degs = map(lambda deg: deg + monomial_deg(sub) - 1, eq_degs)
    quad_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, mon_degs))
    return sum(quad_discrepancies)


def _compute_smd(sub, mlist: list):
    return (monomial_deg(sub) - 1) * len(list(filter(lambda m: monomial_divides(sub, m), mlist)))
