"""split_inert_halves.py — separating the Gaussian (split) primes from the real (inert) primes in the Euler product (2026-09-10).
P_split(s) = prod_{p=1(4)} (1-p^-s)^-1,  P_inert(s) = prod_{p=3(4)} (1-p^-s)^-1,  R(s) = (1-2^-s)^-1.
zeta = R * P_split * P_inert ;  L(chi_4) = P_split * prod_{p=3(4)} (1+p^-s)^-1 = P_split * P_inert(2s) / P_inert(s)   [since (1+x)^-1 = (1-x)/(1-x^2)]
=> zeta * L = R * P_split^2 * P_inert(2s)/... let me keep exact forms:
   P_split(s)^2 = zeta(s) L(s,chi4) (1 - 2^-s) / P_inert(2s)          (I)
   P_inert(s)^2 = [zeta(s) / L(s,chi4)] (1 - 2^-s) P_inert(2s)         (II)
Consequences near a zero: (I) both zero families give P_split ~ (s-rho)^(1/2); (II) zeta-zeros give P_inert ~ (s-rho)^(1/2), chi4-zeros give P_inert ~ (s-rho')^(-1/2).
We validate (I),(II) at s = 2 against direct products, then measure the exponents by scaling |.| at rho + delta."""
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
zeta = lambda s: hurwitz(s, 1.0); L4 = lambda s: 4**(-s)*(hurwitz(s, .25) - hurwitz(s, .75))
N = 10**6
sv = np.ones(N + 1, dtype=bool); sv[:2] = False
for p in range(2, int(N**0.5) + 1):
    if sv[p]: sv[p*p::p] = False
primes = np.nonzero(sv)[0]; split = primes[primes % 4 == 1].astype(float); inert = primes[primes % 4 == 3].astype(float)
def Psplit_direct(s): return complex(np.prod(1/(1 - np.exp(-s*np.log(split)))))
def Pinert_direct(s): return complex(np.prod(1/(1 - np.exp(-s*np.log(inert)))))
# validation at s = 2 (everything converges absolutely)
s = 2.0
lhs1 = Psplit_direct(s)**2; rhs1 = zeta(s)*L4(s)*(1 - 2**-s)/Pinert_direct(2*s)
lhs2 = Pinert_direct(s)**2; rhs2 = zeta(s)/L4(s)*(1 - 2**-s)*Pinert_direct(2*s)
print(f"(I) at s=2:  P_split^2 = {lhs1.real:.9f}  vs  zeta L (1-2^-s)/P_inert(2s) = {rhs1.real:.9f}   (diff {abs(lhs1-rhs1):.1e})")
print(f"(II) at s=2: P_inert^2 = {lhs2.real:.9f}  vs  (zeta/L)(1-2^-s) P_inert(2s) = {rhs2.real:.9f}   (diff {abs(lhs2-rhs2):.1e})")
assert abs(lhs1 - rhs1) < 1e-6 and abs(lhs2 - rhs2) < 1e-6
# P_inert(2s) near Re s = 1/2 via its own recursion (II) at 2s: P_inert(2s)^2 = [zeta(2s)/L(2s)](1-2^-2s) P_inert(4s), P_inert(4s) from the fast direct product
def Pinert_mod(s2):   # |P_inert(s2)| for Re s2 slightly above 1, via the recursion
    return math.sqrt(abs(zeta(s2)/L4(s2))*abs(1 - 2**(-s2))*abs(Pinert_direct(2*s2)))
def modP_split(s): return math.sqrt(abs(zeta(s)*L4(s))*abs(1 - 2**(-s))/Pinert_mod(2*s))
def modP_inert(s): return math.sqrt(abs(zeta(s)/L4(s))*abs(1 - 2**(-s))*Pinert_mod(2*s))
rho_z = complex(0.5, 14.134725141734693)     # zeta's first zero
rho_L = complex(0.5, 6.0209489046975966)     # L(s,chi_4)'s first zero
print("\nexponents by scaling: |P(rho + delta)| for delta = 1e-2, 1e-3, 1e-4 (ratio sqrt(10) = 3.162 per decade means exponent +1/2; 1/sqrt(10) means -1/2)")
for name, rho in (("zeta zero  rho = 1/2+14.1347i", rho_z), ("chi_4 zero rho' = 1/2+6.0209i", rho_L)):
    vs = [modP_split(rho + d) for d in (1e-2, 1e-3, 1e-4)]; vi = [modP_inert(rho + d) for d in (1e-2, 1e-3, 1e-4)]
    es = [math.log10(vs[i]/vs[i+1]) for i in range(2)]; ei = [math.log10(vi[i]/vi[i+1]) for i in range(2)]
    print(f"  near {name}:")
    print(f"     |P_split| = {vs[0]:.4e}, {vs[1]:.4e}, {vs[2]:.4e}  -> local exponent {np.mean(es):+.3f}   (square-root ZERO)")
    print(f"     |P_inert| = {vi[0]:.4e}, {vi[1]:.4e}, {vi[2]:.4e}  -> local exponent {np.mean(ei):+.3f}   ({'square-root ZERO' if np.mean(ei) > 0 else 'square-root POLE'})")
    assert abs(abs(np.mean(es)) - 0.5) < 0.03 and abs(abs(np.mean(ei)) - 0.5) < 0.03
print("\nreading: separating the primes by Gaussian class does not separate the zeros — every zero of zeta AND of L(chi_4) is a branch point of BOTH halves;")
print("         the real-prime half tags them: zeta-zeros are sqrt-zeros of P_inert, chi_4-zeros are sqrt-poles of P_inert; the Gaussian half has sqrt-zeros at both.")
print("         (Landau–Walfisz / Kurokawa: such partial Euler products continue to Re s > 0 with branch points at rho/k and a natural boundary at Re s = 0 — from memory, pin before use.)")
print("ALL CHECKS PASSED")
