from qbee.examples import *
from qbee.util import top_priority
from qbee.experiments import make_benchmark_report



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
    make_benchmark_report({"Default": quad_dom,
                           "Square bound": quad_dom_square_term,
                           "C4 bound": quad_dom_C4_term,
                           "Square + C4 bounds": quad_dom_square_C4_term},
                          n_samples=20, use_multiprocessing=True,
                          title="Quadratization with domination and C4 report")
