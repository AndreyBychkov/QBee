"""
Checks of a given set of new variables, monomial or polynomial, and their reduction.

A set of new variables w = (w_1, ..., w_m), each a polynomial in the states x (not in the inputs u), quadratizes the
system x' = f(x, u) if every right-hand side of x' and of w_i' = grad(w_i) . f is a polynomial of degree at most two
in (x, w, u). With monomial new variables this is a check on single monomials; with polynomial ones (for example
w = (a + b*y)*x) terms can cancel, so the check here is linear algebra:

every monomial of degree <= 2 in (x, u) is a product of two of (1, x, u), so a right-hand side P is quadratic in
(x, w, u) iff the part of P of degree >= 3 lies in the span of the parts of degree >= 3 of the products that involve
a new variable (w_i, w_i x_j, w_i u_k, w_i w_j). The span test is Gaussian elimination modulo the prime p = 2^61 - 1,
with no rounding, and symbolic parameters set to random residues mod p. The check is probabilistic (Schwartz-Zippel):
it answers for generic parameter values, and a wrong answer needs a coincidence of probability about (degree x number
of terms) / p, below 1e-12 for systems of hundreds of terms. A result that holds only for special parameter values
(where terms cancel) is not detected.

    >>> from qbee import functions, is_quadratization, reduce_quadratization
    >>> x, y = functions("x, y")
    >>> system = [(x, y ** 3), (y, x ** 3)]
    >>> is_quadratization(system, [x ** 2, y ** 2, x * y])
    True
    >>> is_quadratization(system, [x ** 2, y ** 2])     # (x**2)' = 2 x y**3 needs x*y
    False
    >>> len(reduce_quadratization(system, [x ** 2, y ** 2, x * y, x ** 2 + y ** 2]))   # drops x**2 + y**2
    3
"""
from __future__ import annotations

import random
from typing import Sequence

import sympy as sp

from .util import INDEPENDENT_VARIABLE

PRIME = 2 ** 61 - 1


class _Echelon:
    """Row echelon basis modulo PRIME over integer columns (rows keyed by pivot column, pivot 1)."""

    def __init__(self):
        self.rows = {}

    def reduce(self, v: dict) -> dict:
        v = dict(v)
        while v:
            piv = min(v)
            row = self.rows.get(piv)
            if row is None:
                return v
            a = v[piv]
            for col, val in row.items():
                nv = (v.get(col, 0) - a * val) % PRIME
                if nv:
                    v[col] = nv
                else:
                    v.pop(col, None)
        return v

    def add(self, v: dict):
        v = self.reduce(v)
        if v:
            piv = min(v)
            inv = pow(v[piv], -1, PRIME)
            self.rows[piv] = {k: val * inv % PRIME for k, val in v.items()}


class QuadratizationCheck:
    """
    Exact quadratization checks of subsets of the new variables `new_vars` of `system` (list of (x, x') pairs, as for
    :func:`qbee.quadratize`; functions without an equation are inputs). Products are computed once, so checking many
    subsets (as :meth:`reduce` does) is cheap.
    """

    def __init__(self, system: Sequence[tuple], new_vars: Sequence[sp.Expr], seed: int = 0):
        self.states = [x for x, _ in system]
        rhs = [sp.sympify(f) for _, f in system]
        funcs = set().union(*(f.atoms(sp.Function) for f in rhs))
        undefined = {f for f in funcs if isinstance(f.func, sp.core.function.UndefinedFunction)}   # as qbee.functions
        self.inputs = sorted((f for f in undefined if f not in self.states), key=str)
        self.gens = self.states + self.inputs
        self.new_vars = [sp.sympify(w) for w in new_vars]
        params = sorted(set().union(*(e.free_symbols for e in rhs + self.new_vars))
                        - {INDEPENDENT_VARIABLE} - set(self.states), key=str)
        rng = random.Random(seed)
        self._values = {p: sp.Integer(rng.randrange(1, PRIME)) for p in params}
        for w in self.new_vars:
            bad = [u for u in self.inputs if w.has(u)]
            if bad:
                raise ValueError(f"new variable {w} contains the inputs {bad}")
        poly = {x: self._poly(f) for x, f in zip(self.states, rhs)}
        W = [self._poly(w) for w in self.new_vars]
        self._rate_high = [(x, self._high(poly[x])) for x in self.states]
        self._new_high = []
        for w in self.new_vars:
            dw = sum((self._poly(sp.diff(w, x)) * poly[x] for x in self.states if w.has(x)), self._poly(0))
            self._new_high.append(self._high(dw))
        lin = [self._poly(1), *[self._poly(g) for g in self.gens]]
        self._pairs = {}
        for i, w in enumerate(W):
            self._pairs[(i, None)] = [h for h in (self._high(w * g) for g in lin) if h]
            for j in range(i + 1):
                self._pairs[(i, j)] = [h for h in [self._high(w * W[j])] if h]

    def _poly(self, e):
        e = sp.expand(sp.sympify(e).xreplace(self._values))
        try:
            p = sp.Poly(e, *self.gens)
        except sp.PolynomialError:
            p = None
        if p is None or not p.domain.is_QQ and not p.domain.is_ZZ:
            raise ValueError(f"not a polynomial in the states and inputs with rational coefficients: {e}")
        return p

    @staticmethod
    def _mod(c):
        r = sp.Rational(c)
        return (int(r.p) % PRIME) * pow(int(r.q) % PRIME, -1, PRIME) % PRIME

    def _high(self, p):
        out = {}
        for mon, c in p.terms():
            if sum(mon) >= 3:
                v = self._mod(c)
                if v:
                    out[mon] = v
        return out

    def failures(self, keep: Sequence[int] | None = None) -> list:
        """The equations (a state, or the index of a new variable) whose right-hand side is not quadratic in
        (x, w_keep, u); empty if the new variables w_keep quadratize the system."""
        keep = sorted(range(len(self.new_vars)) if keep is None else keep)
        prods = [h for i in keep for h in self._pairs[(i, None)]] + \
                [h for i in keep for j in keep if j <= i for h in self._pairs[(i, j)]]
        cols = {}
        for h in prods:
            for mon in h:
                cols.setdefault(mon, len(cols))
        basis = _Echelon()
        for h in prods:
            basis.add({cols[m]: v for m, v in h.items()})
        bad = []
        for name, h in self._rate_high + [(i, self._new_high[i]) for i in keep]:
            if h and (any(m not in cols for m in h) or basis.reduce({cols[m]: v for m, v in h.items()})):
                bad.append(name)
        return bad

    def reduce(self) -> list[int]:
        """Indices of an irredundant subset of the new variables that still quadratizes the system: new variables
        are dropped, last first, in passes until no single one can be dropped. Not necessarily of minimal size."""
        if self.failures():
            raise ValueError("the new variables do not quadratize the system")
        keep, dropped = list(range(len(self.new_vars))), True
        while dropped:
            dropped = False
            for i in reversed(keep):
                trial = [k for k in keep if k != i]
                if not self.failures(trial):
                    keep, dropped = trial, True
        return keep


def is_quadratization(system, new_vars, seed: int = 0) -> bool:
    """Whether the new variables (polynomials in the states) quadratize the system exactly. See the module docstring."""
    return not QuadratizationCheck(system, new_vars, seed).failures()


def reduce_quadratization(system, new_vars, seed: int = 0) -> list:
    """An irredundant subset of the new variables that still quadratizes the system (see QuadratizationCheck.reduce)."""
    check = QuadratizationCheck(system, new_vars, seed)
    return [check.new_vars[i] for i in check.reduce()]


def result_new_variables(result, system) -> list:
    """The new variables of a :class:`QuadratizationResult` of `quadratize(system, input_free=True)` as expressions in
    the functions of `system`, for :func:`is_quadratization` and :func:`reduce_quadratization`. Results with inputs in
    new variables, or of polynomialize_and_quadratize, are not supported."""
    by_name = {}
    for x, f in system:
        for g in sp.sympify(f).atoms(sp.Function) | {x}:
            by_name[g.func.__name__] = g
    out = []
    for eq in result.introduced_variables:
        e = sp.sympify(eq.rhs)
        sub = {s: by_name[s.name] for s in e.free_symbols if s.name in by_name}
        sub.update({g: by_name[g.func.__name__] for g in e.atoms(sp.Function) if g.func.__name__ in by_name})
        out.append(e.xreplace(sub))
    return out
