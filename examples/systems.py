import itertools
import sympy as sp
from qbee import *
from functools import reduce
from random import randrange
from operator import add, mul


def generate_x_plus_1(deg):
    R, x = sp.ring(['x'], sp.QQ)
    return [(x + 1) ** deg]


def generate_dense_polynomial(nvars, deg):
    R, *variables = sp.ring([f"x{i}" for i in range(nvars)], sp.QQ)

    def gen_eq_deg(d: int):
        return reduce(add, [randrange(-10, 10) * reduce(mul, p) for p in itertools.product(variables, repeat=d)]) \
               + randrange(-10, 10)

    def gen_eq():
        return reduce(add, [gen_eq_deg(d) for d in range(1, deg + 1)])

    return [gen_eq() for _ in range(nvars)]


def generate_circular(deg):
    R, x, y = sp.ring(['x', 'y'], sp.QQ)
    return [y ** deg, x ** deg]


def generate_hard(deg):
    R, a, b, c = sp.ring(['a', 'b', 'c'], sp.QQ)
    return [c ** deg + a ** 2 * b ** 2 * c ** 3,
            a ** 2,
            b ** 2]


def generate_hill(k):
    R, h, i, t = sp.ring(['h', 'i', 't'], sp.QQ)
    return [k * i ** 2 * t ** (k - 1),
            -k * i ** 2 * t ** (k - 1),
            R.one]


def generate_lifeware_conjecture(n):
    """
    Generate the system from Conjecture 3.1 from https://hal.inria.fr/hal-02900798v2/document
    Also is called Long Monomial in benchmarks.
    """
    variables = sp.ring([f"x{i}" for i in range(n)], sp.QQ)[1:]
    prod_all = 1
    for i in range(n):
        prod_all *= variables[i] ** 2

    system = []
    for i in range(n):
        system.append(variables[(i + 1) % n] ** 2 + prod_all)
    return system


def generate_cubic_cycle(n):
    variables = sp.ring([f"x{i}" for i in range(n)], sp.QQ)[1:]
    system = []
    for i in range(n):
        system.append(variables[(i + 1) % n] ** 3)
    return system


def generate_cubic_bicycle(n):
    variables = sp.ring([f"x{i}" for i in range(n)], sp.QQ)[1:]
    system = []
    for i in range(n):
        system.append(variables[(i + 1) % n] ** 3 + variables[(i - 1) % n] ** 3)
    return system


def generate_selkov(a, b):
    R, x, y = sp.ring(['x', 'y'], sp.QQ)
    return [-x + a * y + x ** 2 * y,
            b - a * y - x ** 2 * y]


def generate_exponential():
    coef_field = sp.FractionField(sp.QQ, ["Apsi", "Atheta", "bpsi", "btheta", "gamma", "B", "D"])
    Apsi, Atheta, bpsi, btheta, gamma, B, D = coef_field.gens
    R, psi, theta, exp, thinv2, u, du = sp.ring(["psi", "theta", "exp", "thinv2", "u", "u'"], coef_field)
    eqs = [
        Apsi * psi + bpsi * u - D * psi * exp,
        Atheta * theta + btheta * u + B * D * psi * exp
    ]
    eqs.append(gamma * eqs[0] * thinv2 * exp)
    eqs.append(eqs[0] * (-2) * thinv2 * psi)
    eqs.append(du)
    eqs.append(du ** 2)
    return eqs


def generate_cardiovascular_respiratory(n_sites=4, drug_conductance=True, gas_exchange=True, symbolic=False):
    """
    A lumped cardiovascular-respiratory model with piecewise-constant inputs (quadratize with input_free=True).

    * circulation: arterial and venous volumes V_a, V_v; terminal flow G (1 - r_n C_n - r_p C_p)(P_a - P_v) (the drug
      factor only with drug_conductance); a pump driven by a phase oscillator (c, s) at heart rate f, with flow
      u (1 + e C_n) K f (P_v - P_0)(1 - c);
    * a reflex state N, driven by arterial pressure (and tissue O2 and CO2); heart rate f = f_0 + k N + g C_n;
    * two drugs with plasma and effect-site concentrations (C_n, C_p) given by dose-rate inputs; a contractility state u
      and a bleed state B, each driven by a rate input; fluid infusion;
    * with gas_exchange: O2 and CO2 in four compartments (pulmonary, arterial, tissue, venous) carried by the terminal
      flow, with dissolved O2 and n_sites-site haemoglobin binding by mass action, CO2 content linear in its partial
      pressure and the saturation; airway O2 and CO2 fractions, with coefficients w1..w4 and the inspired O2 fraction
      as inputs.

    States: 11, plus 4 (n_sites + 2) + 4 with gas_exchange (39 for the defaults). Inputs: 5, plus 5 with gas_exchange.
    Parameters are symbolic with symbolic=True, otherwise distinct placeholder integers (the search depends only on
    which monomials appear, and runs several times faster without symbolic coefficients).
    """
    Va, Vv, N, c, s, B, u = functions("V_a, V_v, N, c, s, B, u")
    Cn_pl, Cn, Cp_pl, Cp = functions("Cn_pl, C_n, Cp_pl, C_p")
    fluid, bleed, dose_n, dose_p, du = functions("fluid, bleed, dose_n, dose_p, du")
    (Ca, Cv, qa0, qv0, G, r_n, r_p, h_p, K, P0, f0, k_f, g_n, e_n, tau_N, gP, Pref, Vn, CLn, tn, Vp, CLp, tp, b,
     two_pi) = parameters("Ca, Cv, qa0, qv0, G, r_n, r_p, h_p, K, P0, f0, k_f, g_n, e_n, tau_N, gP, Pref, "
                          "Vn, CLn, tn, Vp, CLp, tp, b, two_pi")
    Pa = (Va - qa0) / Ca
    Pv = (Vv - (qv0 - h_p * Cp)) / Cv                        # drug p shifts the venous unstressed volume
    vT = G * (1 - r_n * Cn - r_p * Cp) * (Pa - Pv) if drug_conductance else G * (Pa - Pv)
    f = f0 + k_f * N + g_n * Cn                              # heart rate
    vH = u * (1 + e_n * Cn) * K * f * (Pv - P0) * (1 - c)    # pump flow
    reflex = gP * (Pa - Pref)

    system = []
    if gas_exchange:
        FO2_D, FO2_A, FCO2_D, FCO2_A = functions("FO2_D, FO2_A, FCO2_D, FCO2_A")
        FiO2, w1, w2, w3, w4 = functions("FiO2, w1, w2, w3, w4")
        (gO, gC, Oref, Cref, shunt, Osat, sig, M_O, M_C, kO, kC, Patm, c0, sC, hC, FiCO2, Kl,
         q_pul, q_art, q_tis, q_ven) = parameters("gO, gC, Oref, Cref, shunt, Osat, sig, M_O, M_C, kO, kC, Patm, "
                                                  "c0, sC, hC, FiCO2, Kl, q_pul, q_art, q_tis, q_ven")
        kf = [parameters(f"kf{i}") for i in range(1, n_sites + 1)]
        kb = [parameters(f"kb{i}") for i in range(1, n_sites + 1)]
        comps = ("pul", "art", "tis", "ven")
        dis = {k: functions(f"o_{k}") for k in comps}                      # dissolved O2
        H = {k: [functions(f"H{i}_{k}") for i in range(1, n_sites + 1)] for k in comps}
        Z = {k: functions(f"Z_{k}") for k in comps}                        # CO2 content
        vol = dict(pul=q_pul, art=q_art, tis=q_tis, ven=q_ven)

        def sat(k):
            return sum(i * H[k][i - 1] for i in range(1, n_sites + 1)) / n_sites

        def binding(k):
            """Rates of the binding steps H_(i-1) + O2 <-> H_i, with H_0 = 1 - H1 - ... - H_n."""
            prev = [1 - sum(H[k]), *H[k][:-1]]
            return [kf[i] * dis[k] / sig * prev[i] - kb[i] * H[k][i] for i in range(n_sites)]

        def PCO2(k):
            return (Z[k] - c0 - hC * sat(k)) / sC

        inflow = dict(pul=lambda x: x["ven"], art=lambda x: (1 - shunt) * x["pul"] + shunt * x["ven"],
                      tis=lambda x: x["art"], ven=lambda x: x["tis"])
        flow = dict(pul=(1 - shunt) * vT, art=vT, tis=vT, ven=vT)
        uptake_O = dict(pul=kO * (FO2_A * Patm - dis["pul"] / sig), tis=M_O)   # into the blood (M_O < 0)
        uptake_C = dict(pul=kC * (FCO2_A * Patm - PCO2("pul")), tis=M_C)
        for k in comps:
            rb = binding(k)
            system.append((dis[k], flow[k] * (inflow[k](dis) - dis[k]) / vol[k] - Osat / n_sites * sum(rb)
                           + uptake_O.get(k, 0) / vol[k]))
            for i in range(n_sites):
                h_in = inflow[k]({j: H[j][i] for j in comps})
                nxt = rb[i + 1] if i + 1 < n_sites else 0
                system.append((H[k][i], flow[k] * (h_in - H[k][i]) / vol[k] + rb[i] - nxt))
            system.append((Z[k], (flow[k] * (inflow[k](Z) - Z[k]) + uptake_C.get(k, 0)) / vol[k]))
        reflex += gO * (dis["tis"] + Osat * sat("tis") - Oref) + gC * (Z["tis"] - Cref)
        airway = [
            (FO2_D, w1 * (FiO2 - FO2_D) - w2 * (FO2_D - FO2_A)),
            (FO2_A, w3 * (FO2_D - FO2_A) - w4 * Kl * uptake_O["pul"]),
            (FCO2_D, w1 * (FiCO2 - FCO2_D) - w2 * (FCO2_D - FCO2_A)),
            (FCO2_A, w3 * (FCO2_D - FCO2_A) - w4 * Kl * uptake_C["pul"]),
        ]
    system += [
        (Va, vH - vT),
        (Vv, vT + fluid + b * B - vH),
        (N, (-N + reflex) / tau_N),
        (c, -two_pi * f * s),
        (s, two_pi * f * c),
        (B, bleed),
        (u, du),
        (Cn_pl, (dose_n - CLn * Cn_pl) / Vn),
        (Cn, (Cn_pl - Cn) / tn),
        (Cp_pl, (dose_p - CLp * Cp_pl) / Vp),
        (Cp, (Cp_pl - Cp) / tp),
    ]
    if gas_exchange:
        system += airway
    if symbolic:
        return system
    params = sorted(set().union(*(sp.sympify(rhs).free_symbols for _, rhs in system)) - {INDEPENDENT_VARIABLE},
                    key=str)
    values = {p: i + 2 for i, p in enumerate(params)}
    return [(x, sp.sympify(rhs).subs(values)) for x, rhs in system]
