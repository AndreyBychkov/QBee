from qbee import *
from qbee.examples import *
from copy import deepcopy


def temp_quad(system: EquationSystem, inputs: dict):
    poly_system = polynomialize(system)
    poly_system.print()
    print(poly_system.substitution_equations)

    new_var_name = "z_"
    quad, prunings = PolynomialSystem.from_EquationSystem(poly_system, inputs)
    quad_res = BranchAndBound(quad, pruning_funcs=[pruning_by_best_nvars, pruning_by_quadratic_upper_bound,
                                                   pruning_by_squarefree_graphs] + prunings).quadratize()
    print(quad_res.print(new_var_name))
    quad_system = deepcopy(poly_system)
    quad_monoms = [tuple_to_monom(m, quad_res.system.gen_symbols) for m in quad_res.system.introduced_vars]
    for m in quad_monoms:
        new_var = quad_system.variables.create_variable()
        print(f"{new_var} = {m}")
        quad_system.differential_auxiliary_equation_add(new_var, m)
    quad_system.print()


if __name__ == '__main__':
    A = QQ(1, 2)
    E_a = QQ(3, 2)
    R_u = QQ(1, 1)
    R, c1, c2, c3, c4, w1, w2, w3, w4, w5, w6, T, dT, ddT = ring(
        ['c1', 'c2', 'c3', 'c4', 'w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'T', 'dT', 'ddT'], QQ)
    tmp = -A * w1 * w3 * w5
    system = [
        tmp,
        2 * tmp,
        -tmp,
        -2 * tmp,
        QQ(1, 5) * w1 * w2 * tmp,
        -w2 ** 2 * tmp,
        QQ(13, 10) * w3 * w4 * 2 * tmp,
        -w4 ** 2 * 2 * tmp,
        E_a * R_u * w5 * w6 * dT,
        -2 * T * w6 ** 2,
        dT,
        ddT,
        R.zero
    ]
    input_pruning = partial(pruning_by_excluding_variables, excl_vars=list(map(tuple, [dT, ddT])))
    time_pruning = partial(pruning_by_elapsed_time, start_t=time(), max_t=10)
    quad_system = quadratize(system, pruning_functions=(pruning_by_squarefree_graphs,
                                                        pruning_by_quadratic_upper_bound,
                                                        input_pruning,
                                                        time_pruning),
                             new_vars_name='z')
    print(quad_system)