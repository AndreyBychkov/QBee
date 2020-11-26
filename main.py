from qbee.examples import *

if __name__ == '__main__':
    system = generate_circular(3)
    algo = BranchAndBound(system, aeqd_score, [termination_by_square_bound])
    res = algo.quadratize()
    print(res)