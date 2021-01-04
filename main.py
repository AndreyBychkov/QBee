from qbee.examples import *

if __name__ == '__main__':
    system = generate_circular(4)
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    results = algo.get_optimal_quadratizations()
    for res in results:
        print(res)