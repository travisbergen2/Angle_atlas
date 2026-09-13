"""zero_sum_dictionary.py — receipts for phrasing the instruments as sums over zeros
(Travis, 2026-09-11 11:40).  Uses Odlyzko's first 100,000 zeros (zeros1.txt, on disk).
(1) Bergen-Disc rim angle theta_rho = 2 arctan(1/(2 gamma)) and Li's lambda_n = sum_rho 2 sin^2(n theta/2)
(2) the split/inert race mod 4 (a pure sum over L(s,chi_4) zeros): count and leading period
(3) Legendre's fractional-part error E(x) = x prod(1-1/p) - Phi(x, sqrt x): the seam in zero language
(4) growth exponent of the combination sum = sup Re(rho): zeta vs Davenport-Heilbronn (from mobius_side receipts)
"""
import numpy as np, math, json
R = {}
g = np.loadtxt("/agent/workspace/angle_atlas/zeros1.txt"); g = g[np.isfinite(g)]
T = float(g[-1]); tail_inv_g2 = (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T)      # sum_{gamma>T} 1/gamma^2 (density (1/2pi) log(t/2pi))
theta = 2 * np.arctan(1 / (2 * g))                                                     # rim angle of 1 - 1/rho on the unit circle (RH zeros)
R["rim_angle_first_zero_deg"] = round(float(np.degrees(theta[0])), 5)
R["rim_angle_1/gamma_first"] = round(float(1 / g[0]), 6); R["theta_first"] = round(float(theta[0]), 6)
R["rim_span_deg_100k"] = round(float(np.degrees(theta[0])), 3)
lam_known = {1: 0.0230957089661210, 2: 0.0923457352280475, 3: 0.2076389205, 4: 0.3687904795, 5: 0.5755427145, 6: 0.8275660809}  # Keiper/B-L, from memory
lam = {}
for n in range(1, 7):
    partial = float(np.sum(2 * (1 - np.cos(n * theta))))          # pairs rho, rho-bar: 2(1 - cos n theta) = 4 sin^2(n theta/2)
    est = partial + n * n * tail_inv_g2
    lam[n] = {"partial_100k": round(partial, 7), "tail_n^2*sum1/g^2": round(n * n * tail_inv_g2, 7), "lambda_n_est": round(est, 6), "known": lam_known[n], "diff": f"{est - lam_known[n]:.1e}"}
R["lambda_1_closed_form_1+gamma/2-log(4pi)/2"] = round(1 + 0.5772156649015329 / 2 - math.log(4 * math.pi) / 2, 7)
R["Li_from_rim_angles"] = lam
# (2) the mod-4 race
N = 10**6; isp = np.ones(N + 1, dtype=bool); isp[:2] = False
for p in range(2, 1001):
    if isp[p]: isp[p*p::p] = False
primes = np.nonzero(isp)[0]
race = {}
for x in (10**3, 10**4, 10**5, 10**6):
    q = primes[primes <= x]; race[f"1e{int(math.log10(x))}"] = {"inert_3mod4": int(np.sum(q % 4 == 3)), "split_1mod4": int(np.sum(q % 4 == 1)), "inert_minus_split": int(np.sum(q % 4 == 3) - np.sum(q % 4 == 1))}
R["race_mod4"] = race
gamma_chi4 = 6.020948905                                                                 # first zero of L(s, chi_4), computed from scratch 2026-09-10
R["race_leading_period_in_log_x"] = round(2 * math.pi / gamma_chi4, 4); R["race_leading_period_x_ratio"] = round(math.exp(2 * math.pi / gamma_chi4), 3)
R["race_sqrt_term_from_inert_squares_at_1e6"] = round(math.sqrt(1e6) / math.log(1e6), 1)   # ~ sqrt(x)/log x primes' worth of bias
# (3) Legendre fractional-part error
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, 1001):
    if spf[p] == 0:
        idx = np.arange(p * p, N + 1, p); m = spf[idx] == 0; spf[idx[m]] = p
gpf = np.zeros(N + 1, dtype=np.int64); mu = np.ones(N + 1, dtype=np.int64)
for p in primes:
    gpf[p::p] = p; mu[p::p] *= -1
    if p * p <= N: mu[p * p::p * p] = 0
pi = np.cumsum(isp); E = {}
for x in (10**4, 10**5, 10**6):
    y = math.isqrt(x); d = np.arange(1, x + 1); ok = (mu[1:x + 1] != 0) & (gpf[1:x + 1] <= y); d = d[ok]
    frac_err = float(np.sum(mu[d] * ((x / d) - (x // d))))            # sum mu(d){x/d} = x prod(1-1/p) - Phi
    Phi = int(pi[x] - pi[y] + 1); prod = float(np.prod(1 - 1.0 / primes[primes <= y]))
    E[f"1e{int(math.log10(x))}"] = {"sum_mu(d){x/d}": round(frac_err, 1), "x*prod-Phi": round(x * prod - Phi, 1),
                                    "ratio_to_x/log_x": round(frac_err / (x / math.log(x)), 4)}
R["legendre_fractional_error"] = E; R["limit_ratio_2e^-gamma-1"] = round(2 * math.exp(-0.5772156649015329) - 1, 4)
# (4) exponent = sup Re rho (from mobius_side receipts, 2026-09-10)
R["combination_sum_exponent"] = {"zeta_M(x)_fit_1e4_1e6": 0.466, "zeta_sup_Re_rho_verified_zeros": 0.5,
                                 "DH_B(x)_fit_1e4_1e6": 0.794, "DH_offline_zero_Re": 0.808517}
json.dump(R, open("zero_sum_dictionary_receipts.json", "w"), indent=1)
for k, v in R.items(): print(k, v)
assert abs(lam[1]["lambda_n_est"] - R["lambda_1_closed_form_1+gamma/2-log(4pi)/2"]) < 2e-6
print("ALL ASSERTS PASS")
