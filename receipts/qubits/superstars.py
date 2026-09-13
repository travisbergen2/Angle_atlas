#!/usr/bin/env python3
"""superstars.py (2026-09-12) — does removing 2,3,5,7 make the split/inert reading 'work better'?
E1  sieve usage shares of the smallest prime factor by class over n <= 10^6, with / without p <= 7;
    theory: share(p) = (1/p) * prod_{q<p} (1 - 1/q).
E2  prime-state QFT peaks |S(N/q)| / pi(N) at N = 2^16 vs |mu(q)|/phi(q) (Ramanujan-sum main term),
    with and without the primes 2,3,5,7 deleted.
E3  |S(N/4)| vs |pi(N;1,4) - pi(N;3,4)|.
"""
import numpy as np, json, math

def spf_sieve(N):
    spf = np.zeros(N + 1, dtype=np.int64)
    for p in range(2, int(N ** 0.5) + 1):
        if spf[p] == 0:
            idx = np.arange(p * p, N + 1, p); m = spf[idx] == 0; spf[idx[m]] = p
    n = np.arange(N + 1); spf[(spf == 0) & (n >= 2)] = n[(spf == 0) & (n >= 2)]
    return spf

def cls(p):
    if p == 2: return "ram"
    return "split" if p % 4 == 1 else "inert"

R = {}
# ---------------- E1
N = 10 ** 6; spf = spf_sieve(N); primes = np.flatnonzero(spf == np.arange(N + 1))[1:]  # excludes index 0/1
primes = primes[primes >= 2]
def shares(min_p):
    sel = spf[2:] >= min_p
    tot = sel.sum(); out = {"ram": 0, "split": 0, "inert": 0}
    ps, cnt = np.unique(spf[2:][sel], return_counts=True)
    for p, c in zip(ps, cnt): out[cls(int(p))] += c / tot
    return out
s_all, s_no = shares(2), shares(11)
def theory_share(p, P):
    return (1 / p) * math.prod(1 - 1 / q for q in P if q < p)
P_small = [int(p) for p in primes if p < 200]
th_ratio_window = sum(theory_share(p, P_small) for p in P_small if 11 <= p and cls(p) == "inert") / \
                  sum(theory_share(p, P_small) for p in P_small if 11 <= p and cls(p) == "split")
R["E1"] = dict(shares_all=s_all, ratio_inert_over_split_all=s_all["inert"] / s_all["split"],
               shares_p_ge_11=s_no, ratio_inert_over_split_p_ge_11=s_no["inert"] / s_no["split"],
               theory_ratio_window_11_to_199=th_ratio_window,
               theory_share_3=theory_share(3, P_small), theory_share_5=theory_share(5, P_small), theory_share_7=theory_share(7, P_small))
# ---------------- E2 / E3
M = 2 ** 16; is_p = np.zeros(M, bool); is_p[primes[primes < M]] = True
pi = int(is_p.sum())
def peaks(indicator):
    # INSTRUMENT CORRECTION (first run used FFT bins; 1/q is not on a 2^16 grid for q = 3,5,6,7,9,10,11,12 and the
    # nearest bin smears the peak by a sinc factor — e.g. 0.5*sinc(pi/3)=0.41 at q=3). Evaluate S(1/q) directly.
    ps = np.flatnonzero(indicator).astype(float)
    out = {q: float(abs(np.exp(2j * np.pi * ps / q).sum()) / len(ps)) for q in range(2, 13)}
    S = {"2": np.exp(2j * np.pi * ps / 2).sum(), "4": np.exp(2j * np.pi * ps / 4).sum()}
    return out, S
def mu(q):
    f, m, r = {}, q, 2
    while r * r <= m:
        while m % r == 0: f[r] = f.get(r, 0) + 1; m //= r
        r += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return 0 if any(e > 1 for e in f.values()) else (-1) ** len(f)
def phi(q): return sum(1 for a in range(1, q + 1) if math.gcd(a, q) == 1)
pk_all, S_all = peaks(is_p)
ind_no = is_p.copy(); ind_no[[2, 3, 5, 7]] = False
pk_no, _ = peaks(ind_no)
R["E2"] = {str(q): dict(measured=round(pk_all[q], 5), measured_without_2357=round(pk_no[q], 5),
                        ramanujan=round(abs(mu(q)) / phi(q), 5), mu=mu(q)) for q in range(2, 13)}
p1 = int(np.sum(is_p & (np.arange(M) % 4 == 1))); p3 = int(np.sum(is_p & (np.arange(M) % 4 == 3)))
R["E3"] = dict(N=M, pi=pi, pi_1mod4=p1, pi_3mod4=p3, diff=p1 - p3, abs_S_N_over_4=float(abs(S_all["4"])),
               abs_S_N_over_2=float(abs(S_all["2"])), note="S(1/4) = -1 + i(pi1 - pi3) up to sign convention")
# pre-stated expectations, recorded as pass/fail (never moved after the fact)
R["checks"] = {
    "E1_ratio_all_in_1.9_2.3": bool(1.9 < R["E1"]["ratio_inert_over_split_all"] < 2.3),
    "E1_ratio_p_ge_11_in_1.0_1.2": bool(1.0 <= R["E1"]["ratio_inert_over_split_p_ge_11"] <= 1.2),
    "E2_squarefree_peaks_match_mu_over_phi_within_0.02": all(abs(pk_all[q] - abs(mu(q)) / phi(q)) < 0.02 for q in (2, 3, 5, 6, 7, 10, 11)),
    "E2_nonsquarefree_peaks_below_0.02": all(pk_all[q] < 0.02 for q in (4, 8, 9, 12)),
    "E2_deleting_2357_moves_peaks_less_than_4_over_pi": all(abs(pk_all[q] - pk_no[q]) < 4 / pi + 1e-9 for q in range(2, 13)),
    "E3_S_N_over_4_equals_hypot(1, pi1-pi3)": bool(abs(abs(S_all["4"]) - math.hypot(1, p1 - p3)) < 1e-6),
}
e1 = R["E1"]
print("E1 usage shares of smallest prime factor, n<=1e6")
print("  all p   : ram %.4f  inert %.4f  split %.4f  inert/split = %.3f" % (e1["shares_all"]["ram"], e1["shares_all"]["inert"], e1["shares_all"]["split"], e1["ratio_inert_over_split_all"]))
print("  p>=11   : inert %.4f  split %.4f  inert/split = %.3f   (theory window 11..199: %.3f)" % (e1["shares_p_ge_11"]["inert"], e1["shares_p_ge_11"]["split"], e1["ratio_inert_over_split_p_ge_11"], e1["theory_ratio_window_11_to_199"]))
print("  theory shares 3,5,7: %.4f %.4f %.4f" % (e1["theory_share_3"], e1["theory_share_5"], e1["theory_share_7"]))
print("E2 |S(1/q)|/pi(N), N=2^16, pi=%d   (measured / without 2,3,5,7 / Ramanujan |mu|/phi)" % pi)
for q in range(2, 13):
    v = R["E2"][str(q)]; print("  q=%2d: %.4f  %.4f  %.4f   mu=%2d" % (q, v["measured"], v["measured_without_2357"], v["ramanujan"], v["mu"]))
print("E3", {k: R["E3"][k] for k in ("pi_1mod4", "pi_3mod4", "diff", "abs_S_N_over_4", "abs_S_N_over_2")})
print("CHECKS", json.dumps(R["checks"]))
json.dump(R, open("superstars_receipts.json", "w"), indent=1)
