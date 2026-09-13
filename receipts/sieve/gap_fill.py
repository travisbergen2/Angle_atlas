"""gap_fill.py — receipts for: "the gap between primes is filled with all the
combinations of the primes available; adding a prime adds combinations; that
seems computable."  (Travis, 2026-09-10 16:16)

Every number quoted in the reply is computed here.  numpy only.
"""
import numpy as np, math, json, time

t0 = time.time()
N = 10**6

# ---- smallest prime factor sieve -------------------------------------------
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N ** 0.5) + 1):
    if spf[p] == 0:
        idx = np.arange(p * p, N + 1, p)
        m = spf[idx] == 0
        spf[idx[m]] = p
is_prime = spf == 0
is_prime[:2] = False
n_arr = np.arange(N + 1)
spf[is_prime] = n_arr[is_prime]
primes = np.nonzero(is_prime)[0]
assert len(primes) == 78498

# greatest prime factor and Möbius
gpf = np.zeros(N + 1, dtype=np.int64)
mu = np.ones(N + 1, dtype=np.int64)
for p in primes:
    gpf[p::p] = p
    mu[p::p] *= -1
    if p * p <= N:
        mu[p * p::p * p] = 0
mu[0] = 0

R = {}

# ---- A. every composite in (p_k, p_{k+1}) is a combination of primes <= p_k ---
composite = ~is_prime & (n_arr >= 2)
# for each composite n, the prime before it:
prev_prime = np.zeros(N + 1, dtype=np.int64)
prev_prime[primes] = primes
prev_prime = np.maximum.accumulate(prev_prime)
c = n_arr[composite & (n_arr > 2)]
ratio_gpf = gpf[c] / prev_prime[c]
R["A_max_gpf_over_prev_prime"] = float(ratio_gpf.max())      # must be < 1
R["A_max_spf_needed_below_N"] = int(spf[c].max())            # must be <= 997
R["A_spf_le_sqrt"] = bool(np.all(spf[c] ** 2 <= c))
gaps = np.diff(primes)
R["A_max_gap_below_N"] = int(gaps.max()); R["A_max_gap_at"] = int(primes[gaps.argmax()])
R["A_bertrand_max_ratio_next_over_prev"] = float((primes[1:] / primes[:-1]).max())

# ---- B. the new prime's first exclusive kill is its square -------------------
first_excl = {}
for p in primes[primes <= 1000]:
    k = np.nonzero(spf == p)[0]
    first_excl[int(p)] = int(k[k > p][0])
R["B_first_exclusive_kill_is_square_all_p_le_1000"] = all(v == p * p for p, v in first_excl.items())
R["B_examples"] = {p: first_excl[p] for p in (2, 3, 5, 7, 11, 13, 997)}

# ---- C. Legendre's formula: the signed sum over ALL combinations --------------
def legendre(x, y):
    """sum_{d | P(y)} mu(d) floor(x/d)  — only d <= x contribute."""
    d = n_arr[1:x + 1]
    ok = (mu[1:x + 1] != 0) & (gpf[1:x + 1] <= y)
    d = d[ok]
    return int(np.sum(mu[d] * (x // d))), int(len(d))
pi = np.cumsum(is_prime)
for x, y in ((10**4, 100), (10**6, 1000)):
    val, nterms = legendre(x, y)
    truth = int(pi[x] - pi[y] + 1)
    R[f"C_legendre_x{x}"] = {"value": val, "pi(x)-pi(y)+1": truth, "agree": val == truth,
                             "contributing_combinations": nterms,
                             "all_combinations_2^pi(y)": f"2^{int(pi[y])} = {2.0**int(pi[y]):.3e}"}

# ---- D. the naive combination count: Mertens overshoot -----------------------
for x in (10**4, 10**5, 10**6):
    y = int(math.isqrt(x))
    prod = float(np.prod(1 - 1.0 / primes[primes <= y]))
    naive = x * prod
    truth = int(pi[x] - pi[y] + 1)
    R[f"D_naive_x{x}"] = {"x*prod(1-1/p)": round(naive, 1), "survivors": truth,
                          "ratio": round(naive / truth, 4)}
R["D_limit_2e^-gamma"] = round(2 * math.exp(-0.5772156649015329), 6)

# ---- E. each new prime's share of the integers (as smallest factor) ---------
share = {}
cum = 1.0
for p in (2, 3, 5, 7, 11, 13):
    pred = cum / p
    obs = float(np.mean(spf[2:] == p))
    share[p] = {"predicted": round(pred, 5), "observed_to_1e6": round(obs, 5)}
    cum *= (1 - 1.0 / p)
R["E_share"] = share
R["E_covered_by_primes_le_7"] = round(1 - (1/2)*(2/3)*(4/5)*(6/7), 4)
R["E_covered_by_primes_le_1000"] = round(1 - float(np.prod(1 - 1.0 / primes[primes <= 1000])), 4)
R["E_survivors_fraction_1e6"] = round((pi[N] - pi[1000] + 1) / N, 4)

# ---- F. parity: the combinations cannot tell odd-Omega from even-Omega -------
# Omega via spf
Omega = np.zeros(N + 1, dtype=np.int64)
m = n_arr.copy()
for _ in range(20):
    live = m > 1
    if not live.any(): break
    Omega[live] += 1
    m[live] //= spf[m[live]]
rings = {k: int(np.sum(Omega[2:] == k)) for k in range(1, 8)}
R["F_rings_count_Omega"] = rings
R["F_odd_Omega"] = int(np.sum(Omega[2:] % 2 == 1)); R["F_even_Omega"] = int(np.sum(Omega[2:] % 2 == 0))
lam = np.where(Omega % 2 == 0, 1, -1); lam[0] = 0
R["F_Liouville_sum_L(1e6)"] = int(lam[1:].sum())
R["F_odd_minus_even_over_N"] = round((R["F_odd_Omega"] - R["F_even_Omega"]) / N, 5)

# ---- G. the window again, as 'combinations' ---------------------------------
p0 = 999961
win = {}
for n in range(p0 + 1, 999979 + 1):
    f = []; mm = n
    while mm > 1:
        q = int(spf[mm]); f.append(q); mm //= q
    win[n] = f
R["G_window_factorizations"] = {str(k): v for k, v in win.items()}

R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open("gap_fill_receipts.json", "w"), indent=1)
for k, v in R.items():
    print(f"{k}: {v}")

# hard asserts
assert R["A_max_gpf_over_prev_prime"] < 1
assert R["A_spf_le_sqrt"]
assert R["B_first_exclusive_kill_is_square_all_p_le_1000"]
assert R["C_legendre_x10000"]["agree"] and R["C_legendre_x1000000"]["agree"]
print("ALL ASSERTS PASS")
