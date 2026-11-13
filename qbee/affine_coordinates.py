"""
Polynomial new variables by a linear change of coordinates: the affine sub-expressions of a system's right-hand sides
become coordinates, QBee's monomial search runs in those coordinates, and its new variables, mapped back, are
polynomial new variables of the original system.

In models built from physical sub-expressions, products of affine sub-expressions (a rate f = f0 + k N + g C, a
strength u + e C, a factor 1 - c) quadratize with far fewer new variables than monomials do. Adding such an expression
as an extra variable would make the variables linearly dependent, which QBee's pruning rules do not expect; replacing
a state by it is an invertible linear change of coordinates, so the system stays polynomial of the same degree, the
variables stay independent, and the search and its pruning run unchanged.

    >>> from qbee import functions, parameters
    >>> from qbee.affine_coordinates import quadratize_affine
    >>> N, C, c, s, q = functions("N, C, c, s, q")
    >>> f0, k, g = parameters("f0, k, g")
    >>> f = f0 + k * N + g * C                            # a rate, affine in the states
    >>> system = [(N, -N), (C, -C), (c, -f * s), (s, f * c), (q, -f * c * q)]
    >>> res = quadratize_affine(system)
    >>> f in res.coordinates.values()
    True
    >>> {sp.expand(w) for w in res.new_vars} == {sp.expand(f * c), sp.expand(f * s)}
    True

The coordinates are chosen greedily: the affine sub-expressions (no inputs, at least one state, up to a constant
factor) by how often they appear, each kept if it is linearly independent of those kept and of the states left as
they are; it replaces one of its states (the pivot). The new variables are checked with
:func:`qbee.verification.is_quadratization` on the original system and reduced with
:func:`qbee.verification.reduce_quadratization`.
"""
from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass, field
from typing import Sequence

import sympy as sp

from .quadratization import quadratize
from .util import INDEPENDENT_VARIABLE
from .verification import PRIME, QuadratizationCheck, reduce_quadratization, result_new_variables


def _inputs(system):
    states = [x for x, _ in system]
    funcs = set().union(*(sp.sympify(f).atoms(sp.Function) for _, f in system))
    return [f for f in funcs if isinstance(f.func, sp.core.function.UndefinedFunction) and f not in states]


def _linear_part(expr, states):
    """{state: coefficient} and the constant of an expression affine in the states, or None if it is not affine."""
    e = sp.expand(expr)
    try:
        p = sp.Poly(e, *states)
    except sp.PolynomialError:
        return None
    if p.total_degree() > 1:
        return None
    coeffs = {x: p.coeff_monomial(x) for x in states}
    coeffs = {x: c for x, c in coeffs.items() if c != 0}
    return coeffs, p.coeff_monomial(1)


def affine_subexpressions(system) -> list[tuple[sp.Expr, int]]:
    """
    The affine sub-expressions of the right-hand sides, as written (not expanded): every sum that is affine in the
    states, contains a state and no input, and is not a multiple of a single state. Forms equal up to a constant
    factor are one; returned with their number of occurrences, most frequent first.
    """
    states = [x for x, _ in system]
    inputs = set(_inputs(system))
    counts, forms = Counter(), {}
    for _, f in system:
        for node in sp.preorder_traversal(sp.sympify(f)):
            if not isinstance(node, sp.Add) or (inputs and node.has(*inputs)):
                continue
            lp = _linear_part(node, states)
            if lp is None:
                continue
            coeffs, const = lp
            if not coeffs or (len(coeffs) == 1 and const == 0):
                continue
            # normalise: divide by the coefficient of the first state (in system order)
            lead = next(coeffs[x] for x in states if x in coeffs)
            key = sp.simplify(sp.expand(node / lead))
            counts[key] += 1
            # sympy folds a product's sign into a sum (-f*s holds -f): keep the form without a leading minus
            forms.setdefault(key, -node if node.could_extract_minus_sign() else node)
    return [(forms[k], n) for k, n in counts.most_common()]


def _rank(rows, states, values):
    """Rank (mod PRIME) of the coefficient rows ({state: coeff}) with parameters at the given residues."""
    def mod(c):
        r = sp.Rational(sp.sympify(c).xreplace(values))
        return (int(r.p) % PRIME) * pow(int(r.q) % PRIME, -1, PRIME) % PRIME
    mat = [[mod(r.get(x, 0)) for x in states] for r in rows]
    rank, col = 0, 0
    m = [row[:] for row in mat]
    for col in range(len(states)):
        piv = next((i for i in range(rank, len(m)) if m[i][col]), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        inv = pow(m[rank][col], -1, PRIME)
        m[rank] = [v * inv % PRIME for v in m[rank]]
        for i in range(len(m)):
            if i != rank and m[i][col]:
                a = m[i][col]
                m[i] = [(v - a * w) % PRIME for v, w in zip(m[i], m[rank])]
        rank += 1
    return rank


def choose_coordinates(system, candidates=None, max_coordinates=None, seed: int = 0) -> dict:
    """
    {pivot state: affine form} for the change of coordinates: the candidates (default: affine_subexpressions, most
    frequent first) taken in order, each kept if the kept forms and the states not yet replaced stay linearly
    independent; its pivot is the first of its states (in system order) whose replacement keeps them so.
    """
    states = [x for x, _ in system]
    if candidates is None:
        candidates = [e for e, _ in affine_subexpressions(system)]
    params = sorted(set().union(*(sp.sympify(f).free_symbols for _, f in system), *(sp.sympify(c).free_symbols
                                                                                  for c in candidates))
                    - {INDEPENDENT_VARIABLE} - set(states), key=str)
    rng = random.Random(seed)
    values = {p: sp.Integer(rng.randrange(1, PRIME)) for p in params}
    chosen = {}
    for form in candidates:
        if max_coordinates is not None and len(chosen) >= max_coordinates:
            break
        lp = _linear_part(form, states)
        if lp is None or not lp[0]:
            continue
        coeffs = lp[0]
        for pivot in (x for x in states if x in coeffs and x not in chosen):
            rows = [_linear_part(f, states)[0] for f in chosen.values()] + [coeffs]
            rows += [{x: 1} for x in states if x not in chosen and x != pivot]
            if _rank(rows, states, values) == len(states):
                chosen[pivot] = form
                break
    return chosen


@dataclass
class AffineQuadratization:
    """The result of :func:`quadratize_affine`."""
    coordinates: dict            # {new coordinate (function): affine form in the original states}
    replaced: dict               # {original state: new coordinate it is replaced by}
    system: list                 # the system in the new coordinates
    result: object               # QBee's QuadratizationResult in the new coordinates (None if none was found)
    found: list = field(default_factory=list)      # its new variables, mapped back to the original states
    new_vars: list = field(default_factory=list)   # an irredundant subset of `found` (checked on the original)


def transform(system, chosen: dict, prefix: str = "A"):
    """The system in the coordinates where each pivot state x_p is replaced by A_j = chosen[x_p]."""
    from .util import functions
    states = [x for x, _ in system]
    new = {}
    for j, (pivot, form) in enumerate(chosen.items()):
        new[pivot] = functions(f"{prefix}{j}")
    # x_p in terms of the new coordinates: solve A_j = form_j for the pivots
    eqs = [sp.Eq(new[p], form) for p, form in chosen.items()]
    sol = sp.solve(eqs, list(chosen), dict=True)[0]
    rhs = {x: sp.sympify(f) for x, f in system}
    out = []
    for x in states:
        if x in chosen:
            form = chosen[x]
            rate = sum(sp.diff(form, y) * rhs[y] for y in states if form.has(y))
            out.append((new[x], sp.expand(rate.xreplace(sol))))
        else:
            out.append((x, sp.expand(rhs[x].xreplace(sol))))
    coords = {new[p]: form for p, form in chosen.items()}
    return out, coords, new


def quadratize_affine(system, candidates=None, max_coordinates=None, reduce: bool = True, seed: int = 0,
                      **quadratize_kwargs) -> AffineQuadratization:
    """
    Quadratize `system` (as for :func:`qbee.quadratize`, with input_free=True by default) in coordinates where
    affine sub-expressions replace states (choose_coordinates), and map the new variables back: polynomial new
    variables of the original system, checked with :func:`qbee.verification.is_quadratization` and, with `reduce`,
    reduced to an irredundant subset. `quadratize_kwargs` go to :func:`qbee.quadratize` (pruning_functions,
    calc_upper_bound, ...).
    """
    quadratize_kwargs.setdefault("input_free", True)
    chosen = choose_coordinates(system, candidates, max_coordinates, seed)
    new_system, coords, replaced = transform(system, chosen)
    res = quadratize(new_system, **quadratize_kwargs)
    out = AffineQuadratization(coordinates=coords, replaced=replaced, system=new_system, result=res)
    if not res:
        return out
    back = {y: form for y, form in coords.items()}
    out.found = [sp.expand(w.xreplace(back)) for w in result_new_variables(res, new_system)]
    check = QuadratizationCheck(system, out.found, seed)
    if check.failures():
        raise RuntimeError("the mapped-back new variables do not quadratize the original system")
    out.new_vars = [check.new_vars[i] for i in check.reduce()] if reduce else out.found
    return out
