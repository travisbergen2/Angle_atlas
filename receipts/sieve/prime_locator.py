"""prime_locator.py — is there a law locating unknown primes from the angles of composites around known primes? (2026-09-10)
(1) The angles ARE the sieve: on the sopfr dial, the angle of n is set by its factorization; around a prime p the composites' angles
    are determined by p's residues mod small primes (the wheel); the next prime is the first admissible slot no known prime kills.
(2) Wheel law: every prime > 5 lies in 8 of the 30 residue classes mod 30 (Eratosthenes).
(3) Memorylessness: consecutive gaps are nearly uncorrelated; prime counts in short intervals are ~Poisson (variance ~ mean).
(4) The one real local bias: consecutive primes avoid repeating a residue class (Lemke Oliver–Soundararajan 2016, from memory)."""
import numpy as np, math
N = 10**6
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N**0.5) + 1):
    if spf[p] == 0:
        blk = spf[p*p::p]; blk[blk == 0] = p
spf[spf == 0] = np.arange(N + 1)[spf == 0]; spf[0] = spf[1] = 0
primes = np.nonzero(spf[2:] == np.arange(2, N + 1))[0] + 2
def sopfr(n):
    s = 0
    while n > 1: q = int(spf[n]); s += q; n //= q
    return s
def angle(n): return math.degrees(math.atan(sopfr(n)/n))
# (1) window around a known prime
p = int(primes[-3]); nxt = int(primes[-2])
print(f"(1) window around the known prime p = {p} (next prime {nxt}, gap {nxt - p}):")
print(f"    {'n':>7} {'angle':>8}  smallest factor   residue of p that kills it")
for n in range(p + 1, nxt + 1):
    q = int(spf[n]); a = angle(n)
    kill = f"p ≡ {p % q} (mod {q})  -> p+{n-p} ≡ 0" if n != nxt else "— no small factor: the survivor = next prime (45.00°)"
    print(f"    {n:>7} {a:>8.3f}  {q if n != nxt else '-':>15}   {kill}")
print("    every composite angle is fixed by which small prime divides it, i.e. by p's residues: the angles ARE the sieve. The next prime is what the sieve leaves.")
# (2) wheel
r30 = primes[primes > 5] % 30
classes = sorted(set(r30.tolist()))
print(f"\n(2) wheel law: all {len(r30):,} primes in (5, 10^6] lie in the residue classes {classes} mod 30  [T, Eratosthenes]; 8 of 30 slots admissible")
assert classes == [1, 7, 11, 13, 17, 19, 23, 29]
# (3) memorylessness
gaps = np.diff(primes).astype(float)
c = np.corrcoef(gaps[:-1], gaps[1:])[0, 1]
print(f"\n(3) consecutive gap correlation corr(g_k, g_k+1) over {len(gaps):,} gaps below 10^6 = {c:+.4f}  (near 0: the last gap barely predicts the next)")
# Poisson check: counts of primes in intervals of length L around x ~ 5e5..1e6
x0, L, m = 500000, 200, 2000
starts = x0 + np.arange(m) * ((N - x0 - L) // m)
counts = np.array([np.sum((primes >= a) & (primes < a + L)) for a in starts])
mean, var = counts.mean(), counts.var()
print(f"    primes in {m} intervals of length {L} on [5e5, 1e6]: mean {mean:.2f} (expected L/log x ≈ {L/math.log(7.5e5):.2f}), variance {var:.2f}  -> variance/mean = {var/mean:.2f}")
print("    NOTE (correction on display): a first draft labelled this 'Poisson = 1'. At L = 200 ≈ 15 log x the count is SUB-Poisson — Montgomery–Soundararajan (2004): var/mean ≈ 1 − log L/log x ≈ 0.61 here; Gallagher's Poisson limit (1976) is for L = λ log x, λ fixed. The primes are more regular than random at this scale — a second law, still not a locator.")
# (4) Lemke Oliver–Soundararajan bias
for q in (3, 10):
    res = primes[primes > q] % q
    same = np.mean(res[1:] == res[:-1]); k = len(set(res.tolist()))
    print(f"    consecutive primes with the SAME residue mod {q}: {same:.4f} observed vs {1/k:.4f} if memoryless  -> {'avoid repeating' if same < 1/k else 'no avoidance'}  (Lemke Oliver–Soundararajan 2016)")
print("\nreading: the local pattern gives admissibility (deterministic, exact, as expensive as sieving) and one small explained bias — not a location. The location law is global: the explicit formula, primes <-> zeros.")
print("ALL CHECKS PASSED")
