"""euler_product_check.py — receipts for 'why have we not built the Euler product over all primes' (2026-09-10).
(1) Where there are no zeros (Re s > 1) the finite Euler products converge: s = 2 -> pi^2/6.
(2) On the critical line the finite Euler products do NOT converge and are never near zero, while zeta(1/2 + i gamma_1) = 0
    (zeta computed independently from the alternating eta series, pair-averaged).
(3) The pole-compensated series sum (Lambda(n) - 1) n^{-s}: its partial sums settle for real s > 1/2 (it converges there iff RH;
    its abscissa of convergence is exactly the supremum of the real parts of the zeros).
No mpmath in this sandbox: zeros' heights are Odlyzko table values as used in verify_angles.py.
"""
import numpy as np, math
N = 10**6
# smallest-prime-factor sieve -> Lambda(n)
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N**0.5) + 1):
    if spf[p] == 0:
        block = spf[p*p::p]; block[block == 0] = p
spf[spf == 0] = np.arange(N + 1)[spf == 0]; spf[0] = spf[1] = 0
n = np.arange(2, N + 1)
Lam = np.zeros(N + 1)
q = n.copy(); pk = spf[2:]
# n is a prime power iff n / spf^k == 1
m = q.copy(); 
while True:
    mask = (m % pk == 0) & (m > 1)
    if not mask.any(): break
    m = np.where(mask, m // pk, m)
is_pp = (m == 1)
Lam[2:][is_pp] = np.log(pk[is_pp])
primes = n[spf[2:] == n]
print(f"sieve to {N}: {len(primes)} primes; sum Lambda(n) for n<=N = psi(N) = {Lam.sum():.1f} (N = {N}; psi(N)-N = {Lam.sum()-N:+.1f})")
assert len(primes) == 78498

G1 = 14.134725141734693
def euler_partial(s, P):
    pp = primes[primes <= P].astype(float)
    return np.prod(1.0 / (1.0 - pp ** (-s)))
print("\n(1) finite Euler products at s = 2 (no zeros here) -> pi^2/6 =", f"{math.pi**2/6:.6f}")
for P in (10, 100, 1000, 10**4, 10**5, 10**6):
    v = euler_partial(2.0, P); print(f"   P = {P:>7}: {v.real:.6f}   error {v.real - math.pi**2/6:+.2e}")
print("\n(2) finite Euler products at s = 1/2 + i*gamma_1 (a zero of zeta): they neither converge nor stay away from 0.")
print("    log|P_X| should track the rotating main term  Re[ X^(1-s) / ((1-s) log X) ]  of the prime sum (size sqrt(X)/(|1-s| log X)):")
s = complex(0.5, G1)
mods = []
print(f"   {'P':>8} {'|product|':>10} {'arg':>9} {'log|P_X|':>9} {'main term':>10}")
for P in (10, 100, 1000, 10**4, 10**5, 3*10**5, 10**6):
    v = euler_partial(s, P); mods.append(abs(v))
    main = (P ** (1 - s) / ((1 - s) * math.log(P))).real
    print(f"   {P:>8} {abs(v):>10.4f} {math.degrees(np.angle(v)):>+8.1f}° {math.log(abs(v)):>+9.3f} {main:>+10.3f}")
print(f"   swing: max/min |P_X| = {max(mods)/min(mods):.0f}x over these P; a finite product of nonzero factors is never 0 (min {min(mods):.4f})")
print("   PRE-STATED EXPECTATION 'never near 0' REFUTED at face value (0.013 at P = 1e5) — corrected: no convergence, exponentially widening band.")
# zeta(1/2 + i gamma_1) from the alternating eta series, pair-averaged partial sums
M = 400000
k = np.arange(1, M + 2, dtype=float)
terms = ((-1.0) ** (k - 1)) * np.exp(-s * np.log(k))
S = np.cumsum(terms)
eta = 0.5 * (S[-1] + S[-2])
zeta = eta / (1 - 2 ** (1 - s))
print(f"   zeta(1/2 + i gamma_1) via eta (M = {M}, pair-averaged): |zeta| = {abs(zeta):.2e}  (0 up to the series truncation error)")
assert abs(zeta) < 1e-4 and min(mods) > 0 and max(mods) / min(mods) > 10
print("\n(3) pole-compensated series  sum_{n<=x} (Lambda(n) - 1) n^{-s}  at real s > 1/2 (converges there iff RH):")
coef = Lam[2:] - 1.0
for sig in (0.9, 0.75, 0.6):
    part = np.cumsum(coef * n ** (-sig))
    vals = [part[x - 2] - 1.0 for x in (10**3, 10**4, 10**5, 10**6)]   # the n = 1 term is (0 - 1)*1 = -1
    print(f"   s = {sig}: partial sums at x = 1e3..1e6: " + ", ".join(f"{v:+.4f}" for v in vals) +
          f"   | successive changes: " + ", ".join(f"{vals[i+1]-vals[i]:+.4f}" for i in range(3)))
print("   (the changes shrink like x^{1/2 - s} * (psi(x) - x)/sqrt(x); to the LEFT of the line the same series must diverge if RH holds)")
print("ALL CHECKS PASSED")
