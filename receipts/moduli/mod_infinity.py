"""mod_infinity.py — receipts for 'can we work mod infinity - 1' (2026-09-10).
(1) 2-adic: 1 + 2 + 4 + ... + 2^(n-1) = 2^n - 1 == -1 (mod 2^n): the all-ones number 'infinity - 1' IS -1 in Z_2.
(2) Mersenne moduli 2^k - 1: 2 has multiplicative order exactly k (the bits rotate on a cycle of length k).
(3) A world where RH is a THEOREM: an elliptic curve over F_p. Its zeta function Z(T) = (1 - a_p T + p T^2)/((1-T)(1-pT)),
    a_p = p + 1 - #E(F_p); Hasse (1933): |a_p| <= 2 sqrt(p), i.e. the zeros of the numerator lie on |T| = p^(-1/2).
    In the variable s (T = p^(-s)) that is Re s = 1/2 — a critical LINE that is a CIRCLE: periodic with period 2 pi / log p.
    Curve: y^2 = x^3 + x + 1 (nonsingular for p not dividing the discriminant -16(4+27) = -496 = -2^4 * 31)."""
import math, cmath
# (1)
for n in (4, 8, 16, 32):
    s = sum(2**k for k in range(n)); assert s == 2**n - 1 and (s + 1) % 2**n == 0
print("(1) sum_{k<n} 2^k = 2^n - 1 == -1 (mod 2^n) for n = 4, 8, 16, 32: the binary all-ones number is 2-adically -1  [T]")
# (2)
def order2(m):
    o, x = 1, 2 % m
    while x != 1: x = (2*x) % m; o += 1
    return o
print("(2) order of 2 mod 2^k - 1 for k = 3..12:", [order2(2**k - 1) for k in range(3, 13)], " (= k: the bits rotate)  [T]")
# (3)
def count_points(p, a=1, b=1):
    # #E(F_p) = 1 (infinity) + sum over x of (1 + legendre(x^3 + a x + b))
    n = 1
    for x in range(p):
        r = (x*x*x + a*x + b) % p
        if r == 0: n += 1
        else: n += 2 if pow(r, (p-1)//2, p) == 1 else 0
    return n
print(f"(3) y^2 = x^3 + x + 1 over F_p: a_p = p + 1 - #E, Hasse bound |a_p| <= 2 sqrt(p), Frobenius angle theta_p = arccos(a_p / 2 sqrt p)")
print(f"   {'p':>6} {'#E(F_p)':>8} {'a_p':>5} {'2sqrt(p)':>9} {'|a_p|<=2√p':>11} {'theta_p (deg)':>13} {'|zero T|·√p':>11} {'period 2π/log p':>15}")
worst = 0.0
for p in (5, 7, 11, 13, 101, 1009, 10007):
    if p in (2, 31): continue
    N = count_points(p); a = p + 1 - N; bound = 2*math.sqrt(p)
    assert abs(a) <= bound
    roots = [(a + cmath.sqrt(a*a - 4*p))/2, (a - cmath.sqrt(a*a - 4*p))/2]   # reciprocal roots alpha, alpha-bar of 1 - a T + p T^2
    mod = [abs(r) for r in roots]; worst = max(worst, max(abs(m - math.sqrt(p)) for m in mod))
    theta = math.degrees(math.acos(a/bound))
    print(f"   {p:>6} {N:>8} {a:>5} {bound:>9.3f} {'yes':>11} {theta:>13.2f} {mod[0]/math.sqrt(p):>11.6f} {2*math.pi/math.log(p):>15.4f}")
print(f"   max | |alpha| - sqrt(p) | over these p = {worst:.1e}: every zero on the critical circle — RH for the curve, a theorem (Hasse 1933; Weil 1948 for all curves)  [T]")
print("   In s (T = p^-s): zeros at Re s = 1/2, Im s = ±theta_p/log p + 2πk/log p — finitely many per period: the line is a CIRCLE here.")
print("   Over Q there is no period (12.36 trillion distinct zeros, no repetition): 'mod infinity' never arrives.")
print("ALL CHECKS PASSED")
