from __future__ import annotations

import sympy as sp
from .util import *
from .quadratization import *
from .printer import str_qbee

# The graph from Proposition 4.7
DA_GRAPH = {0: {0, 1, 3}, 1: {1, 2}, 2: {2}, 3: {3}}


def quadratize_dimension_agnostic(system: list[(sp.Symbol, sp.Expr)], 
                                  input_der_orders=None,
                                  input_free=False,
                                  print_intermediate=False,
                                  verbose_output=True):
    """
    Find a dimension-agnostic quadratization for a linearly coupled family of ODE system

    :param system: the symbolic representation of the input family, where Dx is the coupling term for the variables x
    :param input_der_orders: maximal oders for the input variables (if defined, overwrites input_free)
    :param input_free: if an input-free quadratization is sought
    :param print_intermediate: if the quadratization of the auxiliary ODE system from Proposition 4.7 should be printed
    :param verbose_output: if set to True, will print an explanation of how the found dimension-agnostic quadratization works
    :return: dimension-agnistic quadratization of the input family
    """

    odesubs = [{} for _ in range(len(DA_GRAPH))]
    all_funcs = reduce(lambda l, r: l | r, map(lambda e: e.atoms(sp.Function), [e[0] for e in system] + [e[1] for e in system]))

    lhs = [e[0] for e in system]
    inputs = [f for f in all_funcs if str_qbee(f)[0] != 'D' and f not in lhs]
    Dfuncs = [Df for Df in all_funcs if str_qbee(Df)[0] == 'D']
    if len(inputs) > 0:
        print(f"Found inputs {', '.join(map(str_qbee, inputs))}, they will be considered space-independent")

    # duplicating the variables
    for f in all_funcs:
        if f not in Dfuncs and f not in inputs:
            for i in range(len(DA_GRAPH)):
                odesubs[i][f] = functions(str_qbee(f) + "_" + str(i), laurent=f.is_laurent)

    # defining the coupling
    for Df in Dfuncs:
        # checking linearity
        for Dg in Dfuncs:
            if any([e[1].diff(Df).diff(Dg) != 0 for e in system]):
                raise Exception("The family is not linear with respect to the coupling")
        # generating the substitution
        f = [F for F in all_funcs if str_qbee(F) == str_qbee(Df)[1:]]
        if len(f) == 0:
            raise Exception(f"Function {str_qbee(Df)} was interpreted as coupling but the corresponding state {str_qbee(Df)[1:]} was not found")
        f = f[0]
        for i in range(len(DA_GRAPH)):
            odesubs[i][Df] = sum([odesubs[j][f] for j in DA_GRAPH[i]])

    # constructing the system from Proposition 4.7
    ode_system = []
    for lhs, rhs in system:
        for i in range(len(DA_GRAPH)):
            ode_system.append((lhs.subs(odesubs[i]), rhs.subs(odesubs[i])))
    if print_intermediate:
        print("Auxiliary ODE system:")
        for lhs, rhs in ode_system:
            print(f"{str_qbee(lhs)}' = {str_qbee(rhs)}")
        print()

    # computing ODE quadratization
    quadr_system = quadratize(ode_system, 
                              generation_strategy=generation_semidiscretized,
                              input_free=input_free,
                              input_der_orders=input_der_orders)
    if print_intermediate:
        print("Quadratization of the auxiliary ODE system:")
        quadr_system.print()

    # extracting dimension-agnostic quadratization
    da_new_vars_uncoupled = []
    da_new_vars_coupled = []
    uncoupled = {sp.symbols(str_qbee(v)) : sp.symbols(str_qbee(k)) for k, v in odesubs[0].items()}
    coupled = {sp.symbols(str_qbee(v)) : sp.symbols(str_qbee(k) + "_tilda") for k, v in odesubs[1].items()}
    to_avoid = {sp.symbols(str_qbee(v)) for v in odesubs[3].values()}
    for eq in quadr_system.introduced_variables:
        new_var = eq.rhs
        if len(new_var.free_symbols.intersection(uncoupled.keys())) > 0 and len(new_var.free_symbols.intersection(to_avoid)) == 0:
            if len(new_var.free_symbols.intersection(coupled.keys())) > 0:
                da_new_vars_coupled.append(new_var.subs(uncoupled | coupled))
            else:
                da_new_vars_uncoupled.append(new_var.subs(uncoupled))

    if verbose_output:
        print("Every ODE system in the family can be quadratized by adding the following variables")
        if len(da_new_vars_uncoupled) > 0:
            print("\t* For each i, we take variables")
            uncoupled_spec = [v.subs({s: sp.symbols(f"{s}_i") for s in uncoupled.values()}) for v in da_new_vars_uncoupled] 
            print("\t" + ", ".join(map(str, uncoupled_spec)))
        if len(da_new_vars_coupled) > 0:
            ith_block = ", ".join([f"{s}_i'" for s in uncoupled.values() if str(s)[0] != 'D'])
            jth_block = ", ".join([f"{s}_j" for s in uncoupled.values() if str(s)[0] != 'D'])
            print(f"\t* For every i and j such that the variables from the j-th block (that is, {jth_block}) appear in the right-hand sides of the i-th block (that is, in {ith_block}), we take")
            coupled_spec = [
                v.subs({s: sp.symbols(f"{s}_i") for s in uncoupled.values()} | {sp.symbols(f"{s}_tilda"): sp.symbols(f"{s}_j") for s in uncoupled.values()})
                for v in da_new_vars_coupled]
            print("\t" + ", ".join(map(str, coupled_spec)))

    return da_new_vars_uncoupled + da_new_vars_coupled
