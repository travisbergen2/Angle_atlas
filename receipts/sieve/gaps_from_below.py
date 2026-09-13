"""gaps_from_below.py — receipts for Travis's 2026-09-11 06:46 question:
"Can we prove that prime multiples like 2, 3, 5 are used to make composites as
frequently as how far the primes are separated on the other end of the line?
The first prime is used for every even number; the last prime will have a gap
equally as large."
Three layers: (1) exact mean-gap conservation over a period; (2) Mertens/PNT:
the 2e^{-gamma} factor; (3) the small primes' forced gaps: Jacobsthal j(P(y))
and the anatomy of the record prime gaps below 1e6.  numpy only.
"""
import numpy as np, math, json, time
t0 = time.time(); R = {}
GAMMA = 0.5772156649015329; EG = math.exp(GAMMA)

N = 10**6; n = np.arange(N + 1)
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N**0.5) + 1):
    if spf[p] == 0:
        idx = np.arange(p*p, N+1, p); m = spf[idx] == 0; spf[idx[m]] = p
isp = spf == 0; isp[:2] = False; spf[isp] = n[isp]; primes = np.nonzero(isp)[0]; pi = np.cumsum(isp)

# ---- (1) exact conservation over a period: mean gap of survivors = 1/prod(1-1/p) = 1/(1 - sum of shares)
lay1 = {}
for y in (3, 5, 7, 11, 13):
    ps = [int(p) for p in primes if p <= y]
    P = int(np.prod(ps)); phi = int(np.prod([p - 1 for p in ps]))
    shares = []; cum = 1.0
    for p in ps: shares.append(cum / p); cum *= (1 - 1/p)
    lay1[f"y={y}"] = {"P": P, "phi": phi, "mean_gap_P/phi": round(P / phi, 6),
                      "1/(1-sum_shares)": round(1 / (1 - sum(shares)), 6),
                      "e^gamma*log_y": round(EG * math.log(y), 4)}
# wheel gap pattern mod 30 (palindromic) and its Jacobsthal maximum
cop = [a for a in range(1, 31) if math.gcd(a, 30) == 1]
gaps30 = [cop[i+1] - cop[i] for i in range(len(cop)-1)] + [cop[0] + 30 - cop[-1]]
lay1["wheel_mod_30_gaps"] = gaps30; lay1["palindromic"] = gaps30[:-1] == gaps30[:-1][::-1]
R["layer1_exact_conservation"] = lay1

# ---- (2) Mertens vs PNT at the seam y = sqrt(x)
lay2 = {}
for x in (10**4, 10**5, 10**6):
    y = int(math.isqrt(x)); Phi = int(pi[x] - pi[y] + 1)
    prod = float(np.prod(1 - 1.0 / primes[primes <= y].astype(float)))
    lay2[f"x={x}"] = {"observed_mean_gap_x/Phi": round(x / Phi, 3),
                      "sieve_mean_gap_1/prod": round(1 / prod, 3),
                      "Mertens_e^gamma*log_sqrt_x": round(EG * math.log(y), 3),
                      "log_x_local_mean_gap": round(math.log(x), 3),
                      "x/pi(x)": round(x / int(pi[x]), 3),
                      "observed/sieve": round((x / Phi) * prod, 4)}
lay2["limit_observed/sieve"] = round(2 / EG, 4)
R["layer2_Mertens_vs_PNT"] = lay2

# ---- (3a) Jacobsthal j(P(y)): the largest gap the primes <= y can force (scan one period)
jac = {}
for y in (2, 3, 5, 7, 11, 13, 17, 19):
    ps = [int(p) for p in primes if p <= y]; P = int(np.prod(ps))
    cov = np.zeros(P, dtype=bool)
    for p in ps: cov[::p] = True                       # residues divisible by some p (0 mod p)
    surv = np.nonzero(~cov)[0]                        # coprime residues in [0,P)
    g = np.diff(np.concatenate([surv, [surv[0] + P]]))
    j = int(g.max())
    # first prime gap >= j
    pg = np.diff(primes); k = int(np.argmax(pg >= j)) if np.any(pg >= j) else None
    jac[f"y={y}"] = {"P": P, "phi": len(surv), "j(P)": j, "mean_gap": round(P / len(surv), 4),
                     "first_prime_gap>=j_at": (int(primes[k]), int(pg[k])) if k is not None else None}
R["layer3a_Jacobsthal"] = jac

# ---- (3b) anatomy of the record prime gaps below 1e6: who makes the gap?
pg = np.diff(primes); rec = []; best = 0
for i, g in enumerate(pg):
    if g > best:
        best = int(g); p = int(primes[i]); comps = np.arange(p + 1, p + g)
        if len(comps) == 0:                      # the gap 2 -> 3 has no composites
            rec.append({"gap": best, "after_prime": p, "composites": 0}); continue
        s = spf[comps]
        rec.append({"gap": best, "after_prime": p, "composites": int(len(comps)),
                    "covered_by_p<=7": round(float(np.mean(s <= 7)), 3),
                    "covered_by_p<=13": round(float(np.mean(s <= 13)), 3),
                    "covered_by_p<=100": round(float(np.mean(s <= 100)), 3),
                    "largest_spf_needed": int(s.max()),
                    "count_spf>100": int(np.sum(s > 100))})
R["layer3b_record_gaps"] = rec
# expected coverage by 2,3,5,7 of a random stretch: 1 - prod(1-1/p)
R["expected_cover_p<=7"] = round(1 - (1/2)*(2/3)*(4/5)*(6/7), 4)
R["expected_cover_p<=13"] = round(1 - float(np.prod(1 - 1.0/np.array([2,3,5,7,11,13.]))), 4)
R["expected_cover_p<=100"] = round(1 - float(np.prod(1 - 1.0/primes[primes <= 100].astype(float))), 4)
# the 114-gap in full
p = 492113; comps = np.arange(p + 1, p + 114); s = spf[comps]
R["gap114_spf_list"] = [int(v) for v in s]
R["gap114_uncovered_by_p<=13"] = [(int(c), int(v)) for c, v in zip(comps, s) if v > 13]

# ---- (3c) usage counts in the window: how often is p used to make a composite <= x?
x = N
R["usage_2"] = int(x // 2 - 1); R["usage_3"] = int(x // 3 - 1)
R["primes_>x/2_never_used"] = int(pi[x] - pi[x // 2])
R["fraction_of_primes_never_used_in_window"] = round(float(pi[x] - pi[x // 2]) / float(pi[x]), 4)
R["largest_prime_below_1e6"] = int(primes[-1]); R["gap_before_it"] = int(primes[-1] - primes[-2])

R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open("gaps_from_below_receipts.json", "w"), indent=1)
for k, v in R.items():
    if k in ("gap114_spf_list",): print(k, v); continue
    if isinstance(v, dict):
        print(k); [print("   ", kk, vv) for kk, vv in v.items()]
    elif isinstance(v, list):
        print(k); [print("   ", r) for r in v]
    else: print(k, v)
# asserts
for y, r in lay1.items():
    if y.startswith("y="): assert abs(r["mean_gap_P/phi"] - r["1/(1-sum_shares)"]) < 1e-9
assert lay1["palindromic"]
assert jac["y=7"]["j(P)"] == 10 and jac["y=13"]["j(P)"] == 22 and jac["y=19"]["j(P)"] == 34
print("ALL ASSERTS PASS")
