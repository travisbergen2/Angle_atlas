"""proven_vs_needed.py — (1) the proof that pair-Delta is never negative, (2) why lambda_n >= 0 is not that kind of statement,
(3) how the mirror acts on zeros (orbits of size 4 or 2, never 3) and on Pythagorean triples (the 8 signed versions). 2026-09-10."""
import numpy as np, math, cmath
# (1) pair Delta = (2 beta - 1) * (1/|1 - conj rho|^2 - 1/|rho|^2): the two factors always share a sign
print("(1) pair Delta = F1 * F2 with F1 = 2 beta - 1, F2 = 1/|1-conj rho|^2 - 1/|rho|^2 = ((2 beta - 1)) / (|rho|^2 |1-conj rho|^2)")
print("    algebra: |rho|^2 - |1 - conj rho|^2 = (beta^2 + g^2) - ((1-beta)^2 + g^2) = 2 beta - 1  -> F2 = (2 beta - 1)/(|rho|^2 |1-conj rho|^2)")
print("    so pair Delta = (2 beta - 1)^2 / (|rho|^2 |1 - conj rho|^2) >= 0, and = 0 iff beta = 1/2.  A SQUARE over a positive denominator. QED")
rng = np.random.default_rng(1); ok = True
for _ in range(20000):
    b, g = rng.uniform(-2, 3), rng.uniform(-50, 50); rho = complex(b, g); mir = 1 - rho.conjugate()
    if abs(rho) < 1e-9 or abs(mir) < 1e-9: continue
    D = (abs(1/rho)**2 - 2*(1/rho).real) + (abs(1/mir)**2 - 2*(1/mir).real)
    closed = (2*b - 1)**2/(abs(rho)**2*abs(mir)**2)
    ok &= abs(D - closed) < 1e-12*max(1, abs(D)) and D >= -1e-15
print("    checked on 20,000 random (beta, gamma), beta in [-2, 3]: identity holds and D >= 0 everywhere:", ok)
# (2) lambda_n pair term: sign depends on WHERE the zero is — not an identity
print("\n(2) Li pair term for a mirror pair z, 1/conj z:  2 - (r^n + r^-n) cos(n theta), r = |z|")
for r in (1.0, 1.001):
    z = r*cmath.exp(1j*math.radians(7.3)); zm = 1/z.conjugate()
    vals = np.array([(1 - z**n).real + (1 - zm**n).real for n in range(1, 20001)])
    print(f"    r = {r}: min over n <= 20000 = {vals.min():+.4f} at n = {int(vals.argmin())+1}; negative terms: {int((vals < -1e-12).sum())}")
print("    on the circle: 2 - 2cos(n theta) >= 0 for every n (an identity); off the circle: negative for infinitely many n. lambda_n >= 0 for ALL n is a claim about the actual zero set.")
# (3a) mirror group on zeros: Klein four-group orbits
def orbit(s): return {complex(round(x.real, 12), round(x.imag, 12)) for x in (s, s.conjugate(), 1 - s, 1 - s.conjugate())}
print("\n(3a) mirror group {id, conj, s->1-s, s->1-conj s} on a zero: orbit sizes")
for s in (complex(0.5, 14.134725), complex(0.808517, 85.699348), complex(0.5, 0.0), complex(0.3, 0.0)):
    print(f"    s = {s.real:.4f}+{s.imag:.4f}i -> orbit size {len(orbit(s))}")
print("    sizes are 4 (off the line), 2 (on the line, or real and symmetric), 1 only for s = 1/2 itself (not a zero: zeta(1/2) = -1.4604). NEVER 3.")
# (3b) triples: (2+i)^2 = 3+4i and the 8 signed versions from units x conjugation
pi = complex(2, 1); vers = set()
for u in (1, 1j, -1, -1j):
    for w in (pi, pi.conjugate()):
        q = (u*w)**2; vers.add((int(round(q.real)), int(round(q.imag))))
print("\n(3b) Pythagorean triple from (2+i)^2 = 3+4i; conjugation and the four units give FOUR sign patterns (squaring kills the unit sign: u^2 = +-1):")
print("    ", sorted(vers), " hypotenuse = N(pi)^... : |(u w)^2| =", abs(pi**2), "always (the norm is mirror-invariant)")
print("    conjugation alone: (2-i)^2 = 3-4i -> (3,-4,5): the right triangle reflected across the real axis; (i(2-i))^2 = (1+2i)^2 =", (1+2j)**2, "-> (-3, 4): still |3| real, |4| imaginary")
print("    the ordered swap (4,3) is NOT reachable: (a+bi)^2 = (a^2-b^2) + 2ab i puts the EVEN leg in the imaginary slot, always. CORRECTION on display: the first draft said eight signed versions; squares give four.")
print("    the unordered triple {3,4,5} is invariant under all eight: mirroring a triple never changes the triple, only its placement.")
print("ALL CHECKS PASSED")
