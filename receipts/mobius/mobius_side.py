"""mobius_side.py — receipts for adjudicating the relayed chain
   "prime combinations -> Möbius inversion -> 1/zeta(s) -> zero geometry"
(2026-09-10 16:41).  numpy only.  Every quoted number in the reply is here.
"""
import numpy as np, math, json, time, glob, os, urllib.request

t0 = time.time()
N = 10**6
n_arr = np.arange(N + 1)

# ---- sieve: primes, Möbius, Omega ------------------------------------------
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N ** 0.5) + 1):
    if spf[p] == 0:
        idx = np.arange(p * p, N + 1, p); m = spf[idx] == 0; spf[idx[m]] = p
is_prime = spf == 0; is_prime[:2] = False
spf[is_prime] = n_arr[is_prime]
primes = np.nonzero(is_prime)[0]
mu = np.ones(N + 1, dtype=np.int64)
for p in primes:
    mu[p::p] *= -1
    if p * p <= N: mu[p * p::p * p] = 0
mu[0] = 0
Omega = np.zeros(N + 1, dtype=np.int64); m = n_arr.copy()
for _ in range(20):
    live = m > 1
    if not live.any(): break
    Omega[live] += 1; m[live] //= spf[m[live]]
lam = np.where(Omega % 2 == 0, 1, -1); lam[0] = 0

R = {}
sq = np.sqrt(np.maximum(n_arr, 1).astype(float))

# ---- A. Mertens function M(x) = sum mu(n), Liouville L(x) -------------------
M = np.cumsum(mu); L = np.cumsum(lam)
R["A_M_at_powers"] = {f"1e{k}": int(M[10**k]) for k in range(1, 7)}
R["A_L_at_powers"] = {f"1e{k}": int(L[10**k]) for k in range(1, 7)}
ratioM = np.abs(M[2:]) / sq[2:]
R["A_max_|M|/sqrt_x_below_1e6"] = round(float(ratioM.max()), 4)
R["A_argmax"] = int(ratioM.argmax() + 2)
R["A_max_|L|/sqrt_x_below_1e6"] = round(float((np.abs(L[2:]) / sq[2:]).max()), 4)
R["A_L_first_positive_x"] = int(np.nonzero(L[2:] > 0)[0][0] + 2) if np.any(L[2:] > 0) else None

# ---- B. Dirichlet inverse: zeta (must return mu) and Davenport–Heilbronn -----
def dirichlet_inverse(a, N):
    """b with sum_{d|n} a(d) b(n/d) = [n==1]; a[1] must be 1.  O(N log N)."""
    b = np.zeros(N + 1, dtype=float); b[1] = 1.0
    for k in range(1, N // 2 + 1):
        idx = np.arange(2 * k, N + 1, k)
        b[idx] -= a[idx // k] * b[k]
    return b
NB = 10**6
a_zeta = np.ones(NB + 1); a_zeta[0] = 0
b_zeta = dirichlet_inverse(a_zeta, NB)
R["B_inverse_of_zeta_equals_mu"] = bool(np.array_equal(b_zeta[1:].astype(np.int64), mu[1:NB + 1]))
kappa = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
per = np.array([0.0, 1.0, kappa, -kappa, -1.0])       # a(n) by n mod 5 (a(5k)=0)
a_dh = per[np.arange(NB + 1) % 5]; a_dh[0] = 0.0
b_dh = dirichlet_inverse(a_dh, NB)
# sanity: convolution identity on a few n
for n in (2, 3, 6, 10, 30, 97, 1000):
    ds = [d for d in range(1, n + 1) if n % d == 0]
    s = sum(a_dh[d] * b_dh[n // d] for d in ds)
    assert abs(s) < 1e-9, (n, s)
Bdh = np.cumsum(b_dh)
R["B_DH_kappa"] = round(kappa, 9)
R["B_DH_b_first_12"] = [round(float(v), 6) for v in b_dh[1:13]]
R["B_DH_|b(n)|_max_below_1e6"] = round(float(np.abs(b_dh[1:]).max()), 3)
R["B_DH_B(x)/sqrt_x_at_powers"] = {f"1e{k}": round(float(Bdh[10**k] / math.sqrt(10**k)), 3) for k in range(1, 7)}
R["B_zeta_M(x)/sqrt_x_at_powers"] = {f"1e{k}": round(float(M[10**k] / math.sqrt(10**k)), 3) for k in range(1, 7)}
runmax = np.maximum.accumulate(np.abs(Bdh[1:]))
xs = np.array([10**4, 10**5, 10**6]); ys = runmax[xs - 1]
slope = np.polyfit(np.log(xs), np.log(ys), 1)[0]
R["B_DH_running_max_|B|"] = {f"1e{k}": round(float(runmax[10**k - 1]), 2) for k in (3, 4, 5, 6)}
R["B_DH_growth_exponent_fit_1e4_1e6"] = round(float(slope), 3)
R["B_mu_running_max_|M|"] = {f"1e{k}": int(np.maximum.accumulate(np.abs(M[1:]))[10**k - 1]) for k in (3, 4, 5, 6)}

# ---- C. finite reciprocal Euler product = signed combination sum -------------
y = 1000
P_y = primes[primes <= y]
prod2 = float(np.prod(1 - 1.0 / P_y.astype(float) ** 2))
gpf = np.zeros(N + 1, dtype=np.int64)
for p in primes: gpf[p::p] = p
d = n_arr[1:]; ok = (mu[1:] != 0) & (gpf[1:] <= y); d = d[ok]
comb2 = float(np.sum(mu[d] / d.astype(float) ** 2))
R["C_prod_p<=1000_(1-p^-2)"] = round(prod2, 10)
R["C_sum_d<=1e6_mu(d)/d^2"] = round(comb2, 10)
R["C_diff"] = f"{prod2 - comb2:.2e}"
R["C_tail_bound_sum_d>1e6_1/d^2"] = f"{1.0 / N:.1e}"
R["C_6/pi^2"] = round(6 / math.pi ** 2, 10)

# ---- D. the period of the combination space ----------------------------------
theta = float(np.sum(np.log(P_y.astype(float))))
R["D_theta(1000)"] = round(theta, 3)
R["D_log10_P(1000)"] = round(theta / math.log(10), 1)
R["D_phi/P"] = round(float(np.prod(1 - 1.0 / P_y.astype(float))), 6)
R["D_window_fraction_1e6_of_period"] = f"1e{6 - theta / math.log(10):.0f}"
for yy in (10, 100, 1000):
    pr = float(np.prod(1 - 1.0 / primes[primes <= yy].astype(float)))
    R[f"D_mertens_check_y{yy}_prod*e^gamma*log_y"] = round(pr * math.exp(0.5772156649015329) * math.log(yy), 5)

# ---- E. sum over zeros 1/|rho|^2 vs 2*lambda_1 (needs a zeros table) --------
files = glob.glob("/agent/workspace/**/zeros1*", recursive=True)
zpath = files[0] if files else None
if zpath is None:
    try:
        zpath = "/agent/workspace/angle_atlas/zeros1.txt"
        urllib.request.urlretrieve("https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1", zpath)
    except Exception as e:
        R["E_zeros_table"] = f"unavailable: {type(e).__name__}"
        zpath = None
if zpath and os.path.exists(zpath):
    g = np.loadtxt(zpath)
    g = g[np.isfinite(g)]
    lam1 = 1 + 0.5772156649015329 / 2 - math.log(4 * math.pi) / 2
    S = float(np.sum(2.0 / (0.25 + g ** 2)))          # both rho and rho-bar
    T = float(g[-1]); tail = (math.log(T / (2 * math.pi)) + 1) / (math.pi * T)
    R["E_zeros_used"] = int(len(g)); R["E_gamma_max"] = round(T, 3)
    R["E_sum_1/|rho|^2_first_zeros"] = round(S, 7)
    R["E_tail_estimate"] = f"{tail:.2e}"
    R["E_sum_plus_tail"] = round(S + tail, 7)
    R["E_2*lambda_1"] = round(2 * lam1, 7)
    R["E_lambda_1_sum_Re(1/rho)"] = round(float(np.sum(2 * 0.5 / (0.25 + g ** 2))) , 7)

R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open("mobius_side_receipts.json", "w"), indent=1)
for k, v in R.items(): print(f"{k}: {v}")
assert R["B_inverse_of_zeta_equals_mu"]
assert abs(prod2 - comb2) < 2e-6
print("ALL ASSERTS PASS")
