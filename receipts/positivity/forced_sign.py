"""forced_sign.py — 'the arithmetic object cannot leave the circle because this quantity has a forced sign' (2026-09-10 13:07).
Delta(z) = N(z) - Tr(z) for z = 1/rho.  Pointwise: Delta = |z|^2 (1 - 2 beta): NOT sign-forced (positive left of the line, negative right).
Mirror-pair sum: Delta(1/rho) + Delta(1/(1 - conj rho)) = (2 beta - 1)(1/|1-conj rho|^2 - 1/|rho|^2) >= 0, = 0 iff beta = 1/2.
The forced sign is geometry of ANY mirror-symmetric zero set. Exhibit: the Davenport-Heilbronn function (mirror-symmetric, real
coefficients, NO Euler product) has zeros off the line. We locate one from scratch and evaluate the pair Delta there: positive.
L(s, chi) via Hurwitz zeta (Euler-Maclaurin), chi the order-4 character mod 5 with chi(2) = i; kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1)."""
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
chi = {1: 1, 2: 1j, 3: -1j, 4: -1}
def Lchi(s, conj=False):
    return 5**(-s) * sum((np.conj(chi[a]) if conj else chi[a]) * hurwitz(s, a/5) for a in (1, 2, 3, 4))
kappa = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
def DH(s): return 0.5*((1 - 1j*kappa)*Lchi(s) + (1 + 1j*kappa)*Lchi(s, conj=True))
print(f"kappa = {kappa:.6f}; Davenport-Heilbronn f(s) = ((1 - i kappa) L(s,chi) + (1 + i kappa) L(s,chi-bar))/2, coefficients 1, {kappa:.4f}, {-kappa:.4f}, -1, 0 (period 5), real  [T]")
# sanity: DH at s = 2 vs the Dirichlet series summed directly
a = {1: 1.0, 2: kappa, 3: -kappa, 4: -1.0, 0: 0.0}
direct = sum(a[n % 5] / n**2 for n in range(1, 200001))
print(f"check f(2): Hurwitz route {DH(2.0).real:.10f}  vs direct series {direct:.10f}")
assert abs(DH(2.0).real - direct) < 1e-7
def newton(f, s0, it=40):
    s = s0
    for _ in range(it):
        h = 1e-6; d = (f(s + h) - f(s - h))/(2*h); s = s - f(s)/d
    return s
def Delta(z): return abs(z)**2 - 2*z.real
# --- an OFF-line zero of D-H (seed: Spira's 0.8085 + 85.6993 i, from memory — recomputed here) ---
rho = newton(DH, complex(0.81, 85.7))
mirror = 1 - rho.conjugate()
print(f"\noff-line zero found: rho = {rho.real:.6f} + {rho.imag:.6f} i,  |f(rho)| = {abs(DH(rho)):.1e}")
print(f"its mirror 1 - conj(rho) = {mirror.real:.6f} + {mirror.imag:.6f} i is also a zero: |f| = {abs(DH(mirror)):.1e}   (functional-equation symmetry, verified)")
assert abs(DH(rho)) < 1e-8 and abs(DH(mirror)) < 1e-8 and abs(rho.real - 0.5) > 0.2
w1, w2 = 1/rho, 1/mirror
print(f"Delta(1/rho) = {Delta(w1):+.4e}  (beta = {rho.real:.4f} > 1/2: INSIDE the circle, negative);  Delta(1/mirror) = {Delta(w2):+.4e}  (beta = {mirror.real:.4f} < 1/2: OUTSIDE, positive)")
print(f"pair sum = {Delta(w1) + Delta(w2):+.4e}  > 0: the FORCED SIGN HOLDS — and the zero is 0.31 off the line. Forced sign did not force the circle.")
assert Delta(w1) + Delta(w2) > 0
# --- an ON-line zero of D-H for contrast ---
ts = np.arange(4.0, 12.0, 0.01); v = [abs(DH(complex(0.5, t))) for t in ts]; t0 = ts[int(np.argmin(v))]
z_on = newton(DH, complex(0.5, t0))
print(f"\nan on-line zero of D-H for contrast: {z_on.real:.6f} + {z_on.imag:.6f} i, |f| = {abs(DH(z_on)):.1e}; pair Delta = {Delta(1/z_on) + Delta(1/(1 - z_on.conjugate())):+.1e}")
# --- the two positivities ---
print("\nTWO positivities, one free and one not:")
print("  free: pair Delta >= 0 for EVERY mirror-symmetric zero set (zeta's, D-H's, any planted set) — carries no location information.")
print("  dependent: Li/Weil positivity lambda_n >= 0 — TRUE iff all zeros on the line; it FAILS for D-H (off-line zeros make some lambda_n negative):")
# contribution of this one D-H quadruple {rho, mirror, conj rho, conj mirror} to lambda_n = sum (1 - (1 - 1/rho)^n), n up to 3e5
zs = np.array([1 - 1/r for r in (rho, mirror, rho.conjugate(), mirror.conjugate())])
ns = np.arange(1, 300001, dtype=float)
contrib = np.real(4 - np.sum(np.exp(np.outer(ns, np.log(zs))), axis=1))
print("  contribution of this one D-H quadruple to lambda_n: n = 1, 5, 10, 20, 40, 80 ->", ", ".join(f"{contrib[n-1]:+.3f}" for n in (1, 5, 10, 20, 40, 80)), "(positive at small n)")
neg = np.nonzero(contrib < 0)[0]
theta = abs(cmath.phase(zs[0])); eps = abs(abs(zs[1]) - 1)
print(f"  |z| off the circle by {eps:.2e}; arg z = {math.degrees(theta):.3f} deg (period 2pi/theta = {2*math.pi/theta:.0f} in n)")
if len(neg):
    n1 = int(neg[0]) + 1
    print(f"  FIRST NEGATIVE contribution at n = {n1} (value {contrib[n1-1]:+.4f}); minimum over n <= 3e5: {contrib.min():+.3e} at n = {int(np.argmin(contrib))+1}")
    print("  -> the dependent positivity (Li) eventually detects the excursion — slowly, because the zero is high and the offset small; the free positivity never does.")
else:
    print("  no negative contribution below n = 3e5 (it must appear eventually: |z|^n grows without bound)")
print("  CORRECTION on display: the first draft of this receipt asserted 'goes negative' next to values that were all positive; replaced by the computed first-negative n.")
print("ALL CHECKS PASSED")
