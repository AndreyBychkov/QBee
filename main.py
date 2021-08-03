from qbee import *
from qbee.examples import *


def temp_quad(system: EquationSystem, inputs: dict):
    pol = polynomialize(system)
    pol.print()
    print(pol.variables.input)
    print(pol.substitution_equations)

    new_var_name = "z_"
    quad, eqs = PolynomialSystem.from_EquationSystem(pol, inputs, return_equations=True)
    list(map(print, eqs[:len(pol.equations)]))
    quad_res = BranchAndBound(quad, pruning_funcs=[pruning_by_best_nvars, pruning_by_squarefree_graphs]).quadratize()
    print(quad_res.print(new_var_name))
    new_R = eqs[0].ring.drop(*[g for g in eqs[0].ring.gens if str(g) in map(str, pol.variables.parameter)])
    quad_monoms = list(map(lambda m: monom2PolyElem(m, new_R.gens), quad_res.system.introduced_vars))
    # quad_monoms = list(map(lambda m: m.set_ring(eqs[0].ring), quad_monoms))
    quad_system = apply_quadratization(eqs, quad_monoms, system.variables.parameter, new_var_name)

    [print(f"{pol.equations[i].lhs} = {e}") for i, e in enumerate(quad_system[:len(pol.equations)])]
    [print(new_var_name + ("{%d}'" % i) + " = " + str(e)) for i, e in enumerate(quad_system[len(eqs):])]


if __name__ == '__main__':
    # c1, c2, c3, c4 = sp.symbols("c1, c2, c3, c4")
    # A, Ea_Ru = sp.symbols("A, Ea_Ru")
    # T = sp.Symbol("T")
    # tmp = -A * sp.exp(-Ea_Ru / T) * c1 ** 0.2 * c2 ** 1.3
    # dc1, dc2, dc3, dc4 = derivatives([c1, c2, c3, c4])
    # system = EquationSystem([
    #     sp.Eq(dc1, tmp.copy()),
    #     sp.Eq(dc2, 2 * tmp.copy()),
    #     sp.Eq(dc3, -tmp.copy()),
    #     sp.Eq(dc4, -2*tmp.copy())
    # ], [A, Ea_Ru], [T])
    # print("Original system:")
    # system.print()
    # print("="*100)
    # temp_quad(system, {T: 2})
    R, c1, c2, c3, c4, w0, w1, w2, w3, T, dT, ddT = sp.ring("c1, c2, c3, c4, w0, w1, w2, w3, T, T', T''", QQ)
    equations = [
        w0 * c1 * w1 * c2 ** 2 * w3,
        w0 * c1 * w1 * c2 ** 2 * w3,
        w0 * c1 * w1 * c2 ** 2 * w3,
        w0 * c1 * w1 * c2 ** 2 * w3,
        w0**2 * w1 * c2**2 * w3,
        w0 * c1 * c1**2 * c2 * w3,
        w2**2 * dT,
        w2**2 * dT * w3,
        dT,
        ddT,
        R(0)
    ]
    quadratize(equations, new_vars_name="z")
