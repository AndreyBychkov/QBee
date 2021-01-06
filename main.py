import pickle
from qbee.examples import *

def x6x4x3():
    R, x = sp.ring('x', sp.QQ)
    system = PolynomialSystem([x ** 6 + x ** 4 + x ** 3])
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    pickle.dump(results, open('../log/quad_systems.pkl', 'wb'))
    for res in results:
        print(res)

def circ4():
    system = generate_circular(4)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    pickle.dump(results, open('../log/quad_systems.pkl', 'wb'))
    for res in results:
        print(res)


if __name__ == '__main__':
    circ4()
