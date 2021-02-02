from qbee.examples import *
from qbee.vizualization import convex_hull_plot_3d


def x6x4x3():
    R, x = sp.ring('x', sp.QQ)
    system = PolynomialSystem([x ** 6 + x ** 4 + x ** 3])
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)


def circ4():
    system = generate_circular(4)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)


def hard3():
    system = generate_hard(3)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)


if __name__ == '__main__':
    R, x, y = ring(['x', 'y'], QQ)
    system = PolynomialSystem([y ** 4, x ** 2])
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    # algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars, termination_by_newton_polyhedron])
    results = algo.get_quadratizations(4)
    for r in results:
        print(r)