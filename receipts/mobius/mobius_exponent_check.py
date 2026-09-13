#!/usr/bin/env python3
"""mobius_exponent_check.py  (2026-09-11)

Adjudicates the outside-Claude catch on the Möbius "running-max exponent 0.466":
  * recomputes M(x) = sum_{n<=x} mu(n) to 10^6 (independent sieve),
  * reports BOTH candidate statistics with their exact definitions:
      (A) least-squares slope of log(running max |M|) vs log x on the grid
          x = 10^4, 10^5, 10^6   -- the grid mobius_side.py uses for the D-H fit,
      (B) pointwise sup of log|M(x)| / log x over 10^3 <= x <= 10^6,
    plus the record staircase, the global max, and sup |M|/sqrt(x).
Second half: certifies the Davenport-Heilbronn off-line zero 0.808517 + 85.699348 i
by the argument principle (winding number on a small circle, 30-digit arithmetic),
answering the relayed "numerical artifact" claim with a receipt rather than a citation.
"""
import json, math, time, sys
import numpy as np

t0 = time.time()
N = 10**6
R = {}

# ---------- Möbius sieve ----------
mu = np.ones(N + 1, dtype=np.int8)
comp = np.zeros(N + 1, dtype=bool)
for p in range(2, N + 1):
    if comp[p]:
        continue
    if p * p <= N:
        comp[p * p::p] = True
        mu[p * p::p * p] = 0
    mu[p::p] *= -1
mu[0] = 0
M = np.concatenate([[0], np.cumsum(mu[1:], dtype=np.int64)])   # M[x] = M(x)
absM = np.abs(M)

R["M_at_10^k"] = {f"1e{k}": int(M[10**k]) for k in range(1, 7)}
runmax = np.maximum.accumulate(absM)
R["running_max_|M|_at_10^k"] = {f"1e{k}": int(runmax[10**k]) for k in range(1, 7)}

def ls_slope(xs):
    xs = np.array(xs, dtype=float)
    return float(np.polyfit(np.log(xs), np.log(runmax[xs.astype(int)]), 1)[0])

R["(A)_LS_slope_log_runmax_vs_log_x"] = {
    "grid_1e4_1e5_1e6": round(ls_slope([10**4, 10**5, 10**6]), 4),
    "grid_1e3..1e6":    round(ls_slope([10**3, 10**4, 10**5, 10**6]), 4),
    "grid_1e2..1e6":    round(ls_slope([10**2, 10**3, 10**4, 10**5, 10**6]), 4),
    "two_point_1e3_1e6": round(math.log(runmax[10**6] / runmax[10**3]) / math.log(1000), 4),
}

# record staircase (x where |M| sets a new record), x >= 1000
rec_x = [x for x in range(1, N + 1) if absM[x] > (absM[x - 1] if x > 1 else -1) and absM[x] == runmax[x] and (x == 1 or runmax[x] > runmax[x - 1])]
rec_big = [(x, int(absM[x])) for x in rec_x if x >= 1000]
R["record_points_x>=1e3"] = rec_big[-8:]
lx = np.log([x for x, _ in rec_big]); ly = np.log([v for _, v in rec_big])
R["(A')_LS_slope_over_record_points_x>=1e3"] = round(float(np.polyfit(lx, ly, 1)[0]), 4)

# (B) pointwise exponent
xs = np.arange(2, N + 1)
with np.errstate(divide="ignore"):
    e = np.where(absM[2:] >= 2, np.log(np.maximum(absM[2:], 1)) / np.log(xs), 0.0)
for lo in (10, 1000):
    sl = slice(lo - 2, None)
    i = int(np.argmax(e[sl])) + lo
    R[f"(B)_sup_log|M|/log_x_over[{lo},1e6]"] = {"value": round(float(e[i - 2]), 4), "at_x": i, "|M(x)|": int(absM[i])}
gmax = int(absM.max()); where = np.flatnonzero(absM == gmax)
R["global_max_|M|"] = {"value": gmax, "first_x": int(where[0]), "last_x": int(where[-1]), "count_x": int(len(where)),
                       "exponent_at_first_x": round(math.log(gmax) / math.log(int(where[0])), 4)}
R["check_outside_claude_points"] = {"|M(300551)|": int(absM[300551]), "e(300551)": round(math.log(absM[300551]) / math.log(300551), 4),
                                    "|M(926265)|": int(absM[926265]), "e(926265)": round(math.log(max(absM[926265],1)) / math.log(926265), 4)}
ratio = absM[1000:] / np.sqrt(np.arange(1000, N + 1))
j = int(np.argmax(ratio)) + 1000
R["sup_|M|/sqrt(x)_on[1e3,1e6]"] = {"value": round(float(ratio.max()), 4), "at_x": j, "|M(x)|": int(absM[j])}
R["timing_sieve_s"] = round(time.time() - t0, 1)

# ---------- Davenport-Heilbronn certification ----------
try:
    import mpmath as mp
    mp.mp.dps = 30
    kappa = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    def f(s):
        return mp.power(5, -s) * (mp.zeta(s, mp.mpf(1) / 5) + kappa * mp.zeta(s, mp.mpf(2) / 5)
                                  - kappa * mp.zeta(s, mp.mpf(3) / 5) - mp.zeta(s, mp.mpf(4) / 5))
    # validity gate: Hurwitz form vs direct Dirichlet series at s = 2 (coefficients 1, k, -k, -1, 0 periodic)
    a = np.tile(np.array([1.0, float(kappa), -float(kappa), -1.0, 0.0]), N // 5)
    direct2 = float(np.sum(a / np.arange(1, N + 1) ** 2.0))
    R["DH_gate_s=2"] = {"hurwitz": float(mp.re(f(2))), "direct_1e6_terms": round(direct2, 9), "kappa": float(kappa)}
    s0 = mp.findroot(f, mp.mpc("0.808517", "85.699348"))
    s1 = mp.findroot(f, mp.mpc("0.191483", "85.699348"))
    def winding(c, r, n=720):
        vals = [f(c + r * mp.expj(2 * mp.pi * k / n)) for k in range(n)]
        tot = mp.mpf(0); mx = mp.mpf(0)
        for k in range(n):
            d = mp.arg(vals[(k + 1) % n] / vals[k]); tot += d; mx = max(mx, abs(d))
        return float(tot / (2 * mp.pi)), float(min(abs(v) for v in vals)), float(mx)
    w0 = winding(s0, mp.mpf("1e-3")); w1 = winding(s1, mp.mpf("1e-3"))
    R["DH_zero"] = {"s0": mp.nstr(s0, 15), "|f(s0)|": mp.nstr(abs(f(s0)), 3),
                    "winding_r=1e-3": w0[0], "min|f|_on_circle": w0[1], "max_arg_step": w0[2]}
    R["DH_mirror"] = {"s1": mp.nstr(s1, 15), "|s1-(1-conj(s0))|": mp.nstr(abs(s1 - (1 - mp.conj(s0))), 3),
                      "winding_r=1e-3": w1[0], "min|f|_on_circle": w1[1], "max_arg_step": w1[2]}
    R["DH_line_check_Re"] = {"Re(s0)": float(mp.re(s0)), "Re(s0)-1/2": float(mp.re(s0) - mp.mpf(1) / 2)}
except ImportError:
    R["DH"] = "mpmath missing"
R["timing_total_s"] = round(time.time() - t0, 1)

# asserts (registered expectations from the source thread)
assert [R["M_at_10^k"][f"1e{k}"] for k in range(1, 7)] == [-1, 1, 2, -23, -48, 212]
assert [R["running_max_|M|_at_10^k"][f"1e{k}"] for k in range(3, 7)] == [12, 43, 132, 368]
json.dump(R, open("/agent/workspace/angle_atlas/mobius_exponent_check.json", "w"), indent=1)
print(json.dumps(R, indent=1))
