import sympy as sp
from typing import Optional, Collection
from termcolor import colored
from functools import partial
from scipy.spatial import Delaunay, ConvexHull
from sympy.polys.orderings import monomial_key
from .heuristics import aeqd_score
from .quadratization import BranchAndBound, PolynomialSystem, termination_by_best_nvars
from .util import monom2str, get_hull_vertices


def show_convex_hull_report(system: PolynomialSystem, max_depth: Optional[int] = None):
    if max_depth is None:
        algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
        max_depth = algo.quadratize().introduced_vars + algo_hull_size(algo)
        print(f"Upper bound of search is {max_depth}")

    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    quads = filter_duplicates(algo.get_quadratizations(max_depth))

    print(
        f"Order of quadratization: {colored('green', 'green')} if monom is in convex hull and {colored('red', 'red')} if it's not.")
    for quad in sorted(quads, key=PolynomialSystem.new_vars_count):
        is_in_hull = points_in_hull(quad, algo.hull)
        colored_quads = color_quad_in_hull(is_in_hull)

        print(f"{quad.new_vars_count()}: ", end="")
        for q in colored_quads:
            print(q, end=', ')
        print()

    print("=" * 50)

    print(
        f"Order of quadratization: {colored('green', 'green')} if monom is vertex of convex hull and {colored('red', 'red')} if it's not.")
    for quad in sorted(quads, key=PolynomialSystem.new_vars_count):
        belong_to_hull = is_hull_vertices(quad, algo.hull)
        colored_quads = color_quad_in_hull(belong_to_hull)

        print(f"{quad.new_vars_count()}: ", end="")
        for q in colored_quads:
            print(q, end=', ')
        print(f"[{len(list(filter(lambda x: x, belong_to_hull.values())))}/{algo_hull_size(algo)}]")


def filter_duplicates(quads: Collection[PolynomialSystem]):
    res = list()
    intr_vars = list()
    for quad in quads:
        if quad.introduced_vars not in intr_vars:
            res.append(quad)
            intr_vars.append(quad.introduced_vars)
    return res


def color_quad_in_hull(is_in_hull: dict):
    res = list()
    for s, pred in is_in_hull.items():
        if pred:
            res.append(colored(s, 'green'))
        else:
            res.append(colored(s, 'red'))
    return res


def points_in_hull(quad: PolynomialSystem, hull: ConvexHull):
    variables = sorted_vars(quad)
    return dict(zip(
        list(map(partial(monom2str, gens=quad.gen_syms), variables)),
        Delaunay(hull.points).find_simplex(variables) >= 0))


def is_hull_vertices(quad: PolynomialSystem, hull: ConvexHull):
    vertices = get_hull_vertices(hull)
    variables = sorted_vars(quad)
    return dict(zip(
        list(map(partial(monom2str, gens=quad.gen_syms), variables)),
        list(map(lambda m: list(m) in vertices, variables))))


def algo_hull_size(algo: BranchAndBound):
    return len(get_hull_vertices(algo.hull))


def sorted_vars(quad: PolynomialSystem):
    return sorted(quad.introduced_vars,
                  key=lambda m: monomial_key('grlex', reversed(quad.gen_syms))(sp.Monomial(m, quad.gen_syms).as_expr()))
