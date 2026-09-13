#!/usr/bin/env python3
"""qutrit_triangle_mod7.py (2026-09-13) — correcting the 'Eisenstein qutrit triangle' and answering
'taking out the superstars won't even out the legs, will it?'

A. Two-leg triangles (qubit): classes mod 4 and mod 3 to X = 10^6.
   D_pi = pi(x;q,1) - pi(x;q,nonres);  D_theta (log p weights);  D_psi (all prime powers, Lambda weights).
   Superstar removal: recompute D_pi without p in {2,3,5,7}.  Prediction: shift by exactly +1.
   Chebyshev bias: D_pi ~ -sqrt(x)/log x  because every odd prime SQUARE is a residue (1 mod 4, 1 mod 3).
B. Three-leg triangle (qutrit): cubic-residue cosets mod 7 — A={1,6} (cubes), B={3,4}, C={2,5}.
   chi(3^k) = omega^k with 3 a primitive root mod 7.  F3 on (pi_A, pi_B, pi_C) gives
   (S, sum chi(p), sum chibar(p))/sqrt3;  Parseval: pi_A^2+pi_B^2+pi_C^2 = (S^2 + 2|sum chi(p)|^2)/3.
   Squares mod 7 = {1,2,4}: one in each coset -> no square-induced bias between cubic legs (prediction).
C. First zeros of L(s, chi) for the cubic character mod 7 (complex chi; zero set not t-symmetric),
   scanned on Re s = 1/2 for |t| < 20, refined by findroot; regression gate: L(s, chi_3) first zero 8.039737.
Fence: classical number theory + numerics. Nothing bears on RH; GRH for these low zeros is checked, not proven.
"""
import json, math, cmath
import numpy as np
import mpmath as mp

X = 10 ** 6
# ---------------- sieve
is_p = np.ones(X + 1, bool); is_p[:2] = False
for i in range(2, int(X ** 0.5) + 1):
    if is_p[i]: is_p[i * i::i] = False
primes = np.flatnonzero(is_p)
logp = np.log(primes.astype(float))
R = {}

def leg_counts(q, res_class, nonres_class, exclude=()):
    sel = ~np.isin(primes, exclude)
    p = primes[sel]; lp = logp[sel]
    r = p % q
    pi_res = int(np.sum(np.isin(r, res_class))); pi_non = int(np.sum(np.isin(r, nonres_class)))
    th_res = float(lp[np.isin(r, res_class)].sum()); th_non = float(lp[np.isin(r, nonres_class)].sum())
    # psi: add prime powers p^k <= X with weight log p
    ps_res, ps_non = th_res, th_non
    for pp, l in zip(p, lp):
        k = 2; v = pp * pp
        while v <= X:
            rr = v % q
            if rr in res_class: ps_res += l
            elif rr in nonres_class: ps_non += l
            k += 1; v *= pp
    return dict(pi_res=pi_res, pi_non=pi_non, D_pi=pi_res - pi_non, D_theta=round(th_res - th_non, 1), D_psi=round(ps_res - ps_non, 1))

for q, res, non in ((4, (1,), (3,)), (3, (1,), (2,))):
    full = leg_counts(q, res, non); nosup = leg_counts(q, res, non, exclude=(2, 3, 5, 7))
    R[f"A_mod{q}"] = dict(all_primes=full, without_2357=nosup, shift_in_D_pi=nosup["D_pi"] - full["D_pi"],
                          sqrt_x_over_log_x=round(math.sqrt(X) / math.log(X), 1))

# ---------------- B: cubic cosets mod 7
omega = cmath.exp(2j * math.pi / 3)
chi7 = {1: 1, 6: 1, 3: omega, 4: omega, 2: omega ** 2, 5: omega ** 2}   # chi(3^k) = omega^k, 3 primitive root
p7 = primes[primes != 7]; r7 = p7 % 7
piA = int(np.sum((r7 == 1) | (r7 == 6))); piB = int(np.sum((r7 == 3) | (r7 == 4))); piC = int(np.sum((r7 == 2) | (r7 == 5)))
S = piA + piB + piC
sum_chi = piA * 1 + piB * omega + piC * omega ** 2
lhs = piA ** 2 + piB ** 2 + piC ** 2; rhs = (S ** 2 + 2 * abs(sum_chi) ** 2) / 3
squares_mod7 = sorted({(a * a) % 7 for a in range(1, 7)})
R["B_cubic_mod7"] = dict(pi_A_cubes=piA, pi_B=piB, pi_C=piC, S=S, legs_minus_mean=[piA - S / 3, piB - S / 3, piC - S / 3],
                         abs_sum_chi=round(abs(sum_chi), 3), parseval_lhs=lhs, parseval_rhs=round(rhs, 6),
                         parseval_holds=abs(lhs - rhs) < 1e-6, squares_mod7=squares_mod7,
                         note="one square residue per coset -> no square bias between cubic legs; cubes p^3 land only in A (x^(1/3) effect)")
mod4 = R["A_mod4"]["all_primes"]; mod3 = R["A_mod3"]["all_primes"]
R["B_vs_qubit_legs"] = dict(mod4_D_pi=mod4["D_pi"], mod3_D_pi=mod3["D_pi"], mod7_cubic_max_leg_gap=int(max(piA, piB, piC) - min(piA, piB, piC)))

# ---------------- C: zeros of the cubic L-function mod 7 (and L(s,chi_3) regression gate)
mp.mp.dps = 20
def L(s, chi_list):
    return mp.dirichlet(s, chi_list)
chi7_list = [0, 1, omega ** 2, omega, omega, omega ** 2, 1]          # index = residue a mod 7
chi7_list = [mp.mpc(z) if z != 0 else 0 for z in chi7_list]
chi3_list = [0, 1, -1]
def zeros_on_line(chi_list, tmin, tmax, step=0.05):
    ts = np.arange(tmin, tmax, step); vals = [abs(L(mp.mpc(0.5, t), chi_list)) for t in ts]
    out, offline = [], []
    for i in range(1, len(ts) - 1):
        if vals[i] < vals[i - 1] and vals[i] < vals[i + 1] and vals[i] < 0.5:
            try:
                z = mp.findroot(lambda s: L(s, chi_list), mp.mpc(0.5, float(ts[i])))
                if abs(L(z, chi_list)) < 1e-12:
                    if abs(mp.re(z) - mp.mpf(1) / 2) < 1e-8:
                        if all(abs(z - w) > 1e-6 for w in out): out.append(z)
                    elif all(abs(z - w) > 1e-6 for w in offline): offline.append(z)   # e.g. the trivial zero s = 0
            except Exception: pass
    zeros_on_line.offline = [mp.nstr(z, 8) for z in offline]
    return sorted(out, key=lambda z: float(mp.im(z)))
z3 = zeros_on_line(chi3_list, 0.5, 20); off3 = list(zeros_on_line.offline)
z7 = zeros_on_line(chi7_list, -20, 20); off7 = list(zeros_on_line.offline)
n_expected = (20 / math.pi) * math.log(7 * 20 / (2 * math.pi * math.e))   # ~ zero count of L(s,chi), q=7, |t|<20
R["C_zeros"] = dict(L_chi3_first_zeros=[mp.nstr(mp.im(z), 10) for z in z3[:4]], L_chi3_first_zero_re=mp.nstr(mp.re(z3[0]), 12) if z3 else None,
                    regression_8_039737=bool(z3 and abs(float(mp.im(z3[0])) - 8.039737156) < 1e-5),
                    L_cubic7_zeros_t_in_pm20=[mp.nstr(mp.im(z), 10) for z in z7], L_cubic7_max_abs_re_minus_half=mp.nstr(max(abs(mp.re(z) - mp.mpf(1) / 2) for z in z7), 3) if z7 else None,
                    L_cubic7_n_zeros_found=len(z7), L_cubic7_n_expected_riemann_vonMangoldt=round(n_expected, 1),
                    offline_roots_found=dict(chi3=off3, cubic7=off7),
                    note="L(s,chibar) zeros are the conjugates (t -> -t); Dedekind zeta_K (conductor 7 cubic field) = zeta * L(chi) * L(chibar)")
# gate: chi7 is a character (multiplicative on units)
for a in range(1, 7):
    for b in range(1, 7):
        assert abs(chi7[a] * chi7[b] - chi7[(a * b) % 7]) < 1e-12
assert R["B_cubic_mod7"]["parseval_holds"]
for q in (4, 3):
    a = R[f"A_mod{q}"]; f, w = a["all_primes"], a["without_2357"]
    print(f"A mod {q}: pi_res {f['pi_res']} pi_non {f['pi_non']}  D_pi {f['D_pi']}  D_theta {f['D_theta']}  D_psi {f['D_psi']}   | without 2,3,5,7: D_pi {w['D_pi']} (shift {a['shift_in_D_pi']})   sqrt(x)/log x = {a['sqrt_x_over_log_x']}")
b = R["B_cubic_mod7"]
print(f"B cubic mod 7: legs A(cubes) {b['pi_A_cubes']}  B {b['pi_B']}  C {b['pi_C']}  (S {b['S']}); legs-mean {[round(v,1) for v in b['legs_minus_mean']]}; |sum chi(p)| {b['abs_sum_chi']}; Parseval {b['parseval_lhs']} = {b['parseval_rhs']} -> {b['parseval_holds']}; squares mod 7 {b['squares_mod7']}")
print("B vs qubit legs:", R["B_vs_qubit_legs"])
c = R["C_zeros"]
print("C L(chi_3) zeros:", c["L_chi3_first_zeros"], "Re", c["L_chi3_first_zero_re"], "regression 8.039737:", c["regression_8_039737"])
print("C L(cubic mod 7) zeros on Re=1/2, |t|<20:", c["L_cubic7_zeros_t_in_pm20"])
print("   max|Re-1/2|", c["L_cubic7_max_abs_re_minus_half"], " found", c["L_cubic7_n_zeros_found"], " expected ~", c["L_cubic7_n_expected_riemann_vonMangoldt"], " off-line roots:", c["offline_roots_found"])
json.dump(R, open("qutrit_triangle_mod7_receipts.json", "w"), indent=1, default=str)
