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
    c1, c2, c3, c4 = sp.symbols("c1, c2, c3, c4")
    A, Ea, Ru = sp.symbols("A, Ea, Ru")
    T = sp.Symbol("T")
    tmp = -A * sp.exp(-Ea / (Ru * T)) * c1 ** 0.2 * c2 ** 1.3
    dc1, dc2, dc3, dc4 = derivatives([c1, c2, c3, c4])
    system = EquationSystem([
        sp.Eq(dc1, tmp.copy()),
        sp.Eq(dc2, 2 * tmp.copy()),
        sp.Eq(dc3, -tmp.copy()),
        sp.Eq(dc4, -2 * tmp.copy())
    ], [A, Ea, Ru], [T])
    print("Original system:")
    system.print()
    print("=" * 100)
    temp_quad(system, {T: 2})
