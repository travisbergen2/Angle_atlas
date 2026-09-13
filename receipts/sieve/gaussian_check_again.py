"""gaussian_check_again.py — "Now what if we separate the Gaussian primes and
check again" (Travis, 2026-09-11 07:43).  Redo the usage-from-below / gaps-
above layers with the rational primes separated into 2 (ramified), split
(p = 1 mod 4, two Gaussian primes of norm p) and inert (p = 3 mod 4, one
Gaussian prime of norm p^2); then run the sieve INSIDE Z[i] by norm.

PRE-STATED before running (falsifiable):
  P1  the seam ratio observed/sieve in Z[i] at X = 1e6 lies within 0.02 of
      the rational value 1.0336 (both tend to 2e^{-gamma} = 1.1229; the field
      constant 4/pi cancels).
  P2  inert primes carry a larger usage share than split primes below 1000
      (3 is inert).
  P3  consecutive primes in the same class mod 4 have a LARGER mean gap than
      class-switching pairs, and the same-class fraction is < 1/2
      (4 | gap  <=>  same class; gap 2 always switches).
numpy only.
"""
import numpy as np, math, json, time
t0 = time.time(); R = {"prestated": ["P1 seam ratio Z[i] within 0.02 of 1.0336", "P2 inert usage share > split", "P3 same-class mean gap > switching mean gap; same-class fraction < 1/2"]}
GAMMA = 0.5772156649015329; EG = math.exp(GAMMA)

# ---- rational sieve --------------------------------------------------------
N = 10**6; n = np.arange(N + 1)
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N**0.5) + 1):
    if spf[p] == 0:
        idx = np.arange(p*p, N+1, p); m = spf[idx] == 0; spf[idx[m]] = p
isp = spf == 0; isp[:2] = False; spf[isp] = n[isp]; primes = np.nonzero(isp)[0]; pi = np.cumsum(isp)
cls = lambda p: "ram" if p == 2 else ("split" if p % 4 == 1 else "inert")

# ---- A. usage shares by class (smallest prime factor) ------------------------
s_cls = {"ram": 0.0, "split": 0.0, "inert": 0.0}; cum = 1.0
for p in primes[primes <= 1000]:
    s_cls[cls(int(p))] += cum / p; cum *= (1 - 1.0 / p)
obs = {c: 0.0 for c in s_cls}
sp = spf[2:]; obs["ram"] = float(np.mean(sp == 2))
obs["split"] = float(np.mean((sp % 4 == 1) & (sp <= 1000))); obs["inert"] = float(np.mean((sp % 4 == 3) & (sp <= 1000)))
R["A_usage_share_p<=1000"] = {c: {"predicted": round(s_cls[c], 5), "observed_1e6": round(obs[c], 5)} for c in s_cls}
R["A_counts_p<=1000"] = {"split": int(np.sum(primes[primes <= 1000] % 4 == 1)), "inert": int(np.sum(primes[primes <= 1000] % 4 == 3))}

# ---- B. consecutive primes: class transitions mod 4 and gaps ----------------
odd = primes[primes > 2]; g = np.diff(odd); same = (odd[1:] % 4) == (odd[:-1] % 4)
assert np.array_equal(same, g % 4 == 0)              # 4 | gap <=> same class
R["B_same_class_fraction"] = round(float(np.mean(same)), 4)
R["B_mean_gap_same_class"] = round(float(g[same].mean()), 3)
R["B_mean_gap_switch"] = round(float(g[~same].mean()), 3)
R["B_ratio_same/switch"] = round(float(g[same].mean() / g[~same].mean()), 3)
trans = {}
for a in (1, 3):
    for b in (1, 3):
        sel = (odd[:-1] % 4 == a); trans[f"{a}->{b}"] = round(float(np.mean(odd[1:][sel] % 4 == b)), 4)
R["B_transition_matrix_mod4"] = trans

# ---- C. record gaps WITHIN each class, and the 114-gap killers by class ------
for a, name in ((1, "split"), (3, "inert")):
    q = primes[primes % 4 == a]; gq = np.diff(q); k = int(gq.argmax())
    R[f"C_{name}_class_below_1e6"] = {"count": int(len(q)), "mean_gap": round(float(gq.mean()), 3), "2*log(1e6)": round(2 * math.log(1e6), 3),
                                      "record_gap": int(gq[k]), "after": int(q[k])}
comps = np.arange(492114, 492113 + 114); s = spf[comps]
R["C_gap114_killers_by_class"] = {"ram(2)": int(np.sum(s == 2)), "split": int(np.sum(s % 4 == 1)), "inert": int(np.sum(s % 4 == 3))}
R["C_gap114_killers_by_class_excluding_2"] = {"split": int(np.sum((s % 4 == 1))), "inert": int(np.sum(s % 4 == 3)), "note": "of the 56 composites not killed by 2"}

# ---- D. Jacobsthal of class primorials --------------------------------------
def jac(ps):
    P = int(np.prod(ps)); cov = np.zeros(P, dtype=bool)
    for p in ps: cov[::p] = True
    surv = np.nonzero(~cov)[0]; gg = np.diff(np.concatenate([surv, [surv[0] + P]]))
    return P, len(surv), int(gg.max())
split_ps = [int(p) for p in primes if p % 4 == 1][:5]; inert_ps = [int(p) for p in primes if p % 4 == 3][:6]
R["D_jacobsthal_split_first_k"] = {}; R["D_jacobsthal_inert_first_k"] = {}
for k in range(2, 6):
    P, phi, j = jac(split_ps[:k]); R["D_jacobsthal_split_first_k"][k] = {"primes": split_ps[:k], "P": P, "mean_gap": round(P / phi, 4), "j": j}
for k in range(2, 7):
    P, phi, j = jac(inert_ps[:k]); R["D_jacobsthal_inert_first_k"][k] = {"primes": inert_ps[:k], "P": P, "mean_gap": round(P / phi, 4), "j": j}

# ---- E. the sieve INSIDE Z[i]: lattice points of norm <= X, Gaussian primes of norm <= y
X = 10**6; y = 1000
B = int(math.isqrt(X)); a, b = np.meshgrid(np.arange(-B, B + 1), np.arange(-B, B + 1), indexing="ij")
a = a.ravel(); b = b.ravel(); nm = a * a + b * b; keep = (nm <= X) & (nm > 0); a = a[keep]; b = b[keep]; nm = nm[keep]
npts = int(len(a)); R["E_lattice_points_0<N<=X"] = npts; R["E_pi*X"] = round(math.pi * X)
killed = np.zeros(npts, dtype=bool); first_norm = np.zeros(npts, dtype=np.int64)
# Gaussian primes of norm <= y, in norm order
gp = [(2, "ram", 2)]                                  # (norm, type, p)
for p in primes[primes <= y]:
    p = int(p)
    if p % 4 == 1: gp.append((p, "split", p))
    elif p % 4 == 3 and p * p <= y: gp.append((p * p, "inert", p))
gp.sort()
prod = 1.0
for norm, typ, p in gp:
    if typ == "ram":
        k = ((a + b) % 2 == 0); prod *= (1 - 1 / 2)
    elif typ == "inert":
        k = (a % p == 0) & (b % p == 0); prod *= (1 - 1 / (p * p))
    else:
        r = next(r for r in range(1, p) if (r * r + 1) % p == 0)
        k = (((a + b * r) % p) == 0) | (((a - b * r) % p) == 0); prod *= (1 - 1 / p) ** 2
    new = k & ~killed; first_norm[new] = norm; killed |= k
surv = int(np.sum(~killed))
# independent count: 4 * (#prime ideals with y < norm <= X) + 4 units
split_big = int(np.sum((primes % 4 == 1) & (primes > y) & (primes <= X)))
inert_big = int(np.sum((primes % 4 == 3) & (primes > math.isqrt(y)) & (primes * primes <= X)))
ideal_count = 2 * split_big + inert_big
R["E_survivors_observed"] = surv; R["E_4*(prime_ideals_y<N<=X)+4units"] = 4 * ideal_count + 4
R["E_gaussian_Mertens_prod_N<=1000"] = round(prod, 6)
R["E_Rosen_(4/pi)e^-gamma/log_y"] = round((4 / math.pi) * math.exp(-GAMMA) / math.log(y), 6)
R["E_sieve_prediction_pts*prod"] = round(npts * prod)
R["E_seam_ratio_observed/sieve"] = round(surv / (npts * prod), 4)
# CONVENTION CORRECTION (kept on display): P1 as written quoted the rational
# value 1.0336, which is the sieve/observed COUNT ratio (= observed/sieve mean-GAP
# ratio).  This block computes observed/sieve COUNTS, whose rational twin at the
# same X is 78,331/80,965 = 0.9675.  Both conventions are reported; the
# pre-statement is judged twice: as literally written, and in matched convention.
R["E_rational_seam_count_ratio_same_X_observed/sieve"] = round(78331 / 80965, 4)
R["E_rational_seam_ratio_as_quoted_in_P1_(sieve/observed)"] = 1.0336
R["E_gaussian_sieve/observed"] = round((npts * prod) / surv, 4)
R["E_limit_sieve/observed"] = round(2 / EG, 4); R["E_limit_observed/sieve"] = round(EG / 2, 4)
# Gaussian Mertens convergence (products only need rational primes)
def gm(yy):
    pr = 1.0
    for p in primes[primes <= yy]:
        p = int(p)
        if p == 2: pr *= 0.5
        elif p % 4 == 1: pr *= (1 - 1 / p) ** 2
        elif p * p <= yy: pr *= (1 - 1 / (p * p))
    return pr
R["E_gaussian_Mertens_check_prod*log_y*(pi/4)e^gamma"] = {yy: round(gm(yy) * math.log(yy) * (math.pi / 4) * EG, 5) for yy in (10**3, 10**4, 10**5, 10**6)}
# Gaussian usage shares by first killer norm
pred = {2: 0.5, 5: 0.5 * (1 - (4 / 5) ** 2), 9: 0.5 * (4 / 5) ** 2 * (1 / 9), 13: 0.5 * (4 / 5) ** 2 * (8 / 9) * (1 - (12 / 13) ** 2)}
R["E_gaussian_shares"] = {nn: {"predicted": round(v, 5), "observed": round(float(np.mean(first_norm == nn)), 5)} for nn, v in pred.items()}

# ---- verdicts on the pre-statements --------------------------------------------
R["verdict_P1_as_written"] = abs(R["E_seam_ratio_observed/sieve"] - 1.0336) < 0.02       # False: convention error in the pre-statement (mine)
R["verdict_P1_matched_convention"] = abs(R["E_seam_ratio_observed/sieve"] - 78331 / 80965) < 0.02
R["P1_difference_matched"] = round(abs(R["E_seam_ratio_observed/sieve"] - 78331 / 80965), 4)
R["verdict_P2"] = s_cls["inert"] > s_cls["split"]
R["verdict_P3"] = (R["B_mean_gap_same_class"] > R["B_mean_gap_switch"]) and (R["B_same_class_fraction"] < 0.5)
R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open("gaussian_check_again_receipts.json", "w"), indent=1, default=str)
for k, v in R.items(): print(k, v)
assert surv == 4 * ideal_count + 4
print("ALL ASSERTS PASS")
