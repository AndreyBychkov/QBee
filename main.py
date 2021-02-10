from qbee.examples import *


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
    system = generate_lifeware_conjecture(6)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    algo.inside_newton_polyhedral_upper_bound()
    quad_res = algo.quadratize()
    print(quad_res)
