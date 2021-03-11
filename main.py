from qbee.examples import *
from qbee.util import top_priority, calc_Lie_derivative, generalized_variables_dict
from qbee.experiments import make_benchmark_report, get_examples


def quad(system):
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    start_t = time()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_square_term(system):
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars, termination_by_square_bound])
    start_t = time()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_C4_term(system):
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars, termination_by_C4_bound])
    start_t = time()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_square_C4_term(system):
    algo = BranchAndBound(system, aeqd_score,
                          [termination_by_best_nvars, termination_by_C4_bound, termination_by_square_bound])
    start_t = time()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_dom(system):
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars])
    start_t = time()
    algo.domination_upper_bound()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_dom_square_term(system):
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars, termination_by_square_bound])
    start_t = time()
    algo.domination_upper_bound()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_dom_C4_term(system):
    algo = BranchAndBound(system, aeqd_score, [termination_by_best_nvars, termination_by_C4_bound])
    start_t = time()
    algo.domination_upper_bound()
    algo.quadratize()
    end_t = time()
    return end_t - start_t


def quad_dom_square_C4_term(system):
    algo = BranchAndBound(system, aeqd_score,
                          [termination_by_best_nvars, termination_by_C4_bound, termination_by_square_bound])
    start_t = time()
    algo.domination_upper_bound()
    algo.quadratize()
    end_t = time()
    return end_t - start_t

if __name__ == '__main__':
    R, x = sp.ring('x', sp.QQ)
    system = [x**5,]
    res = quadratize(system)
    print(res)