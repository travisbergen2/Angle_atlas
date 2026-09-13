"""product_positivity.py — positive forms derivable from the Euler product WITHOUT RH, and how far they reach (2026-09-10).
(A) 3-4-1: 3 + 4 cos(phi) + cos(2 phi) = 2 (1 + cos phi)^2 >= 0, combined with Lambda(n) >= 0 (the product's coefficient positivity),
    gives 3 log zeta(s) + 4 Re log zeta(s+it) + Re log zeta(s+2it) >= 0 for s > 1, i.e. zeta(s)^3 |zeta(s+it)|^4 |zeta(s+2it)| >= 1.
    Consequence: no zero on Re s = 1 (PNT), and with derivative bounds a zero-free region s > 1 - 1/(R log t), R = 5.573412 (explicit, 2018).
(B) A Gram (Hardy-space) form: K(s,t) = zeta(s + conj t) = sum_n n^-s conj(n^-t) is positive semi-definite for Re s, Re t > 1/2 — RH-free — and blind to zeros.
(C) The Weil form's positivity (Li's lambda_n >= 0) is the target; NOT derivable here — it is RH.
L-values via Hurwitz zeta (Euler–Maclaurin)."""
import numpy as np, math, cmath
from math import factorial
B2k = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330]
def hurwitz(s, a, Nt=1200, m=10):
    n = np.arange(Nt, dtype=float); S = np.sum(np.exp(-s*np.log(n + a))); x = Nt + a
    S += x**(1 - s)/(s - 1) + 0.5*x**(-s); rf = s
    for k in range(1, m + 1):
        if k > 1: rf *= (s + 2*k - 3)*(s + 2*k - 2)
        S += B2k[k-1]/factorial(2*k) * rf * x**(-s - 2*k + 1)
    return S
zeta = lambda s: hurwitz(s, 1.0)
# (A) the square
phis = np.linspace(0, 2*math.pi, 100001); tri = 3 + 4*np.cos(phis) + np.cos(2*phis)
assert np.all(tri >= -1e-12) and np.allclose(tri, 2*(1 + np.cos(phis))**2)
print("(A) 3 + 4cos(phi) + cos(2phi) = 2(1 + cos phi)^2 >= 0 on [0, 2pi]  [T]  — a SQUARE turns coefficient positivity into an inequality")
print("    the Euler product gives log zeta(s) = sum_n Lambda(n)/(log n) n^-s with EVERY coefficient >= 0 (Re s > 1)  [T]")
print("    => zeta(s)^3 |zeta(s+it)|^4 |zeta(s+2it)| >= 1 for all real s > 1, t:   (Mertens 1898; de la Vallee Poussin 1899)")
print(f"    {'sigma':>7} {'t':>8} {'zeta(s)^3 |zeta(s+it)|^4 |zeta(s+2it)|':>40}")
worst = np.inf
for sig in (1.5, 1.1, 1.02, 1.005, 1.001):
    for t in (14.134725, 21.022040, 30.0, 85.699348):
        v = abs(zeta(sig))**3 * abs(zeta(complex(sig, t)))**4 * abs(zeta(complex(sig, 2*t)))
        worst = min(worst, v)
        if t in (14.134725, 85.699348): print(f"    {sig:>7} {t:>8.3f} {v:>40.6f}")
print(f"    minimum over the grid = {worst:.6f}  (>= 1 always)  [T ✓]")
assert worst >= 1 - 1e-9
print("    mechanism: if zeta(1 + i t0) = 0 then near sigma -> 1+: zeta(sigma)^3 ~ (sigma-1)^-3 but |zeta(sigma+it0)|^4 ~ (sigma-1)^4, product -> 0 < 1. Contradiction. No zeros on Re s = 1 -> PNT.")
R = 5.573412   # explicit zero-free region constant, live-sourced today (arXiv 1809.03134 abstract): no zeros with sigma >= 1 - 1/(R log t), t >= 3
print(f"\n    how far the product's positivity reaches (explicit region, R = {R}): sigma > 1 - 1/(R log t)")
for t in (1e2, 1e6, 1e12, 3e12):
    print(f"       t = {t:.0e}: zero-free for sigma > {1 - 1/(R*math.log(t)):.5f}   (RH would need 1/2; the gap is {1 - 1/(R*math.log(t)) - 0.5:.5f})")
print("    the method lives where the product converges (Re s > 1) and leaks inward by ~1/log t; Vinogradov–Korobov (1958) reach 1 - c/(log t)^(2/3)(log log t)^(1/3) with exponential sums; nothing reaches 1/2 + epsilon.")
# (B) Gram form
pts = [complex(0.6, 0.0), complex(0.7, 5.0), complex(0.8, 14.134725), complex(1.2, -3.0), complex(0.55, 30.0)]
K = np.array([[zeta(a + b.conjugate()) for b in pts] for a in pts])
K = 0.5*(K + K.conj().T)
ev = np.linalg.eigvalsh(K)
print(f"\n(B) Gram form K_ij = zeta(s_i + conj s_j) on five points with Re s > 1/2: eigenvalues min {ev.min():.4e}, max {ev.max():.4e} -> positive definite  [T]")
assert ev.min() > 0
print("    RH-free positivity (a Gram matrix of the vectors n^-s), but it lives in Re s > 1/2 and cannot see the zeros: it is blind by domain, not by weakness.")
# (C)
print("\n(C) the form whose positivity IS RH: Weil's explicit-formula form / Li's lambda_n / Paper 14's A(m). Built from the product side; positivity unproven; failing D–H (rule 10) but an equivalence, not a theorem.")
print("ALL CHECKS PASSED")
