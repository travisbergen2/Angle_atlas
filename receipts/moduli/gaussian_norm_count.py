"""gaussian_norm_count.py — is the counting variable x a norm? (2026-09-10)
In Z the size of n is |n|; in Z[i] the size of z is the norm z*zbar (a product of conjugates, the shape of (+i)(-i) = 1).
The Gaussian prime number theorem counts Gaussian primes by NORM <= x: pi_{Q(i)}(x) ~ li(x), and under the Riemann
hypothesis for zeta_{Q(i)}(s) = zeta(s) L(s, chi_4) the error is O(sqrt(x) log x) — von Koch's shape with x a norm.
Receipt: count Gaussian primes of norm <= x exactly from a rational-prime sieve and compare with li(x).
Gaussian primes: 1+i (norm 2, one prime up to units); each split p = 1 mod 4 gives two primes (norm p); each inert
p = 3 mod 4 stays prime with norm p^2 (so counts iff p^2 <= x)."""
import numpy as np, math
N = 10**6
sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
for p in range(2, int(N**0.5) + 1):
    if sieve[p]: sieve[p*p::p] = False
primes = np.nonzero(sieve)[0]
def li(x):  # li(x) = Ei(ln x) = gamma + ln t + sum t^k/(k k!)
    t = math.log(x); s = 0.5772156649015329 + math.log(t); term = 1.0
    for k in range(1, 200):
        term *= t / k; s += term / k
        if term / k < 1e-16: break
    return s
print(f"{'x':>8} {'pi_Q(i)(x)':>11} {'li(x)':>10} {'diff':>8} {'sqrt(x)log x':>13} | {'pi(x)':>7} {'pi-li':>8}")
for e in (3, 4, 5, 6):
    x = 10**e
    p_le = primes[primes <= x]
    split = int(np.sum(p_le % 4 == 1)); inert_sq = int(np.sum((primes[primes <= int(math.isqrt(x))] % 4) == 3))
    piK = 2 * split + inert_sq + 1                     # + the ramified prime 1+i
    print(f"{x:>8} {piK:>11} {li(x):>10.1f} {piK - li(x):>+8.1f} {math.sqrt(x)*math.log(x):>13.0f} | {len(p_le):>7} {len(p_le)-li(x):>+8.1f}")
# the anchor: (+i)(-i) = 1 is where every counting function starts
print("\n(+i)(-i) =", complex(0,1)*complex(0,-1), "= N(i); pi(1) = psi(1) = 0; li(1) = -inf; |x^rho| = sqrt(x) = 1 at x = 1.")
print("Gaussian primes by norm: 1+i (norm 2), split p -> two primes of norm p, inert p -> one prime of norm p^2.")
print("Check at 10^6: split =", int(np.sum(primes % 4 == 1)), "inert p <= 1000 =", int(np.sum((primes[primes <= 1000] % 4) == 3)), "-> pi_Q(i)(10^6) = 2*split + inert + 1")
