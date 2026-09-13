"""eisenstein.py — the Gaussian pictures rebuilt in Z[w], w = e^{2 pi i/3} (2026-09-10).
Units: six (the 6th roots of unity), rotation by 60 degrees; conjugation folds; fundamental angle domain [0, 30] degrees.
Ramified prime 1 - w (norm 3): 3 = -w^2 (1-w)^2; arg(1 - w) = -30 deg == 30 deg mod 60: it sits ON the fold edge.
Norm N(a + b w) = a^2 - a b + b^2; trace = 2a - b. p is a norm iff p = 3 or p = 1 mod 3 (chi_3 = (p/3)); p = 2 mod 3 inert.
Trace = norm  <=>  2 Re z = |z|^2  <=>  |z - 1| = 1: the SAME circle as in Z[i]; its lattice points are 1 + (units): six here.
Dedekind zeta zeta_{Q(w)}(s) = zeta(s) L(s, chi_3) = Euler product over Eisenstein primes by norm.
L-functions via Hurwitz zeta (Euler-Maclaurin), validated on zeta's first zero and L(s, chi_4)'s first zero."""
import numpy as np, math, cmath
from math import factorial
w = cmath.exp(2j*math.pi/3)
deg = math.degrees
# ---------- units, fold, ramified prime ----------
units = [(-w)**k for k in range(6)]
unit_angles = sorted(round(deg(cmath.phase(u)) % 360, 6) for u in units)
assert unit_angles == [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
rp = 1 - w
assert abs(abs(rp)**2 - 3) < 1e-12 and abs((1 - w)**2 - (-3*w)) < 1e-12 and abs(-w**2*(1-w)**2 - 3) < 1e-12
print("units of Z[w]: six, at", unit_angles, "degrees; fundamental angle domain after rotation + conjugation: [0, 30] degrees  [T]")
print(f"ramified prime 1 - w: norm {abs(rp)**2:.0f}, arg {deg(cmath.phase(rp)):+.1f} deg == {deg(cmath.phase(rp)) % 60:.1f} deg mod 60: ON the fold edge; 3 = -w^2 (1-w)^2  [T]")
# ---------- norm rule ----------
N = 10**6
sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
for p in range(2, int(N**0.5) + 1):
    if sieve[p]: sieve[p*p::p] = False
primes = np.nonzero(sieve)[0]
def rep(p):  # p = x^2 + 3 y^2  ->  pi = x + y sqrt(-3) = (x + y) + 2 y w
    for y in range(0, math.isqrt(p // 3) + 1):
        x2 = p - 3*y*y; x = math.isqrt(x2)
        if x*x == x2: return x, y
    return None
small = [int(p) for p in primes if p < 100]
sol = {p: rep(p) for p in small}
assert all((sol[p] is not None) == (p == 3 or p % 3 == 1) for p in small)
print("norm rule checked on primes < 100: p is a norm iff p = 3 or p = 1 mod 3  [T, Fermat/Euler]")
print("  split (norms):", [p for p in small if p % 3 == 1], " inert (not norms):", [p for p in small if p % 3 == 2])
print("  e.g. 7 = N(3 + 2w) [x=2,y=1], 13 = N(", (rep(13)[0] + rep(13)[1], 2*rep(13)[1]), "), 3 = N(1 - w)")
# ---------- trace = norm circle: lattice points ----------
pts = []
for a in range(-4, 5):
    for b in range(-4, 5):
        z = a + b*w
        if abs(abs(z - 1) - 1) < 1e-9: pts.append((a, b, z))
assert len(pts) == 6
print("trace = norm circle |z - 1| = 1 in Z[w]: SIX lattice points (1 + each unit):", [f"{a}+{b}w" for a, b, _ in pts], " [T]")
print("  their reciprocals sit on Re s = 1/2 at heights", sorted(round(-(1/z).imag, 4) for a, b, z in pts if abs(z) > 0), "and 0 -> infinity;  Z[i] had four points (0, 2, 1±i) at heights 0, ±1/2")
# ---------- Hurwitz zeta (Euler-Maclaurin) and the L-functions ----------
B2k = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330]
def hurwitz(s, a, Nt=1000, m=10):
    n = np.arange(Nt, dtype=float)
    S = np.sum(np.exp(-s*np.log(n + a)))
    x = Nt + a
    S += x**(1 - s)/(s - 1) + 0.5*x**(-s)
    rf = s
    for k in range(1, m + 1):
        if k > 1: rf *= (s + 2*k - 3)*(s + 2*k - 2)
        S += B2k[k-1]/factorial(2*k) * rf * x**(-s - 2*k + 1)
    return S
def zeta(s): return hurwitz(s, 1.0)
def L4(s): return 4**(-s)*(hurwitz(s, 0.25) - hurwitz(s, 0.75))
def L3(s): return 3**(-s)*(hurwitz(s, 1/3) - hurwitz(s, 2/3))
def first_zero(f, lo, hi, step=0.002):
    ts = np.arange(lo, hi, step); v = [abs(f(complex(0.5, t))) for t in ts]
    i = int(np.argmin(v)); a, b = ts[max(i-1, 0)], ts[min(i+1, len(ts)-1)]
    for _ in range(60):
        c = b - 0.6180339887*(b - a); d = a + 0.6180339887*(b - a)
        if abs(f(complex(0.5, c))) < abs(f(complex(0.5, d))): b = d
        else: a = c
    t = 0.5*(a + b); return t, abs(f(complex(0.5, t)))
tz, vz = first_zero(zeta, 10.0, 16.0); t4, v4 = first_zero(L4, 4.0, 8.0); t3, v3 = first_zero(L3, 4.0, 12.0)
assert abs(tz - 14.134725141734693) < 1e-6 and abs(t4 - 6.0209489046975966) < 1e-6
print(f"validation: zeta first zero {tz:.9f} (|zeta| {vz:.1e}); L(s,chi_4) first zero {t4:.9f} (|L| {v4:.1e}) — both match the record  [T ✓]")
print(f"L(s, chi_3) first zero on the line, from scratch: t = {t3:.9f}  (|L| = {v3:.1e});  at that height |zeta| = {abs(zeta(complex(0.5, t3))):.4f}, |L(chi_4)| = {abs(L4(complex(0.5, t3))):.4f}")
# ---------- Dedekind zeta at s = 2 ----------
L3_2_direct = float(np.sum(1/(3*np.arange(0, 2_000_000, dtype=float) + 1)**2 - 1/(3*np.arange(0, 2_000_000, dtype=float) + 2)**2))
L3_2 = L3(2.0).real
assert abs(L3_2 - L3_2_direct) < 1e-9
target = (math.pi**2/6)*L3_2
print(f"L(2, chi_3) = {L3_2:.12f} (direct sum {L3_2_direct:.12f}); zeta(2) L(2, chi_3) = {target:.10f}")
for X in (10**3, 10**4, 10**5, 10**6):
    pp = primes[primes <= X]; split = pp[pp % 3 == 1].astype(float); inert = pp[(pp % 3 == 2) & (pp*pp <= X)].astype(float)
    prod = (1/(1 - 3.0**-2)) * np.prod((1 - split**-2.0)**-2) * np.prod((1 - inert**-4.0)**-1)
    print(f"  Eisenstein Euler product over norms <= {X:>7}: {prod:.10f}   error {prod - target:+.2e}")
assert abs(prod - target) < 1e-6
# ---------- angle census of split primes ----------
split_p = [int(p) for p in primes if p % 3 == 1]
folded = []; cos6 = {}
for p in split_p:
    x, y = rep(p); th = deg(math.atan2(y*math.sqrt(3), x)) % 60.0
    if th > 30.0: th = 60.0 - th
    folded.append(th)
    if p in (7, 13, 19, 31, 37, 43):
        pi_ = complex(x, y*math.sqrt(3)); c6 = (pi_**6).real/p**3
        cos6[p] = (round(c6*p**3), p**3, c6)
folded = np.array(folded); fs = np.sort(folded)/30.0
ks = float(np.max(np.abs(fs - (np.arange(1, len(fs)+1)/len(fs)))))
print(f"split primes <= 10^6: {len(split_p)} (inert: {int(np.sum(primes % 3 == 2))}); folded angles in [0, 30] deg: Kolmogorov distance to uniform = {ks:.5f} (Hecke equidistribution shadow)  [T]/measured")
print("  cos 6 theta_pi = Re(pi^6)/p^3, exact rationals:", {p: f"{n}/{d}" for p, (n, d, c) in cos6.items()})
np.save("eisenstein_folded.npy", folded)
import json
json.dump({"t3": t3, "t4": t4, "tz": tz, "L3_2": L3_2, "ks": ks, "n_split": len(split_p), "pts": [f"{a}+{b}w" for a, b, _ in pts]}, open("eisenstein_receipts.json", "w"))
print("ALL CHECKS PASSED")
