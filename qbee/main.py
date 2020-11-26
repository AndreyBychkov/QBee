from examples import *

if __name__ == '__main__':
    system = selkov(2, 2)
    algo = BranchAndBound(system, aeqd_score, [termination_by_square_bound])
    res = algo.quadratize()
    print(res)