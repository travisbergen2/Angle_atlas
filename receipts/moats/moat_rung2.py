"""moat_rung2.py — decide the k = 4 (and then √18, √20, √26) Gaussian moats at
larger radius, point-based (no grid arrays), so radius 6000–10000 fits in memory.
Usage: python3 moat_rung2.py RADIUS
Closure rule: the component of 1+i under steps <= k is CLOSED if its farthest
point satisfies |z| <= RADIUS - k (every Gaussian prime in the disk is known).
"""
import sys, numpy as np, math, json, time
from scipy import sparse
from scipy.sparse.csgraph import connected_components
t0 = time.time(); B = int(sys.argv[1]) if len(sys.argv) > 1 else 6000; X = B * B; W = 2 * B + 1
# rational primes to X
isp = np.ones(X + 1, dtype=bool); isp[:2] = False
for p in range(2, math.isqrt(X) + 1):
    if isp[p]: isp[p*p::p] = False
t_sieve = time.time() - t0
# Gaussian primes as points (a, b), norm <= X: a^2+b^2 prime, or inert p on the axes
bb = np.arange(-B, B + 1, dtype=np.int64); A_list = []; B_list = []
for a in range(-B, B + 1):
    nm = a * a + bb * bb; ok = nm <= X; sel = np.zeros_like(ok); sel[ok] = isp[nm[ok]]
    if a == 0: sel |= isp[np.abs(bb)] & (np.abs(bb) % 4 == 3)
    if abs(a) % 4 == 3 and isp[abs(a)]: sel[B] = True            # (a, 0) inert
    b_sel = bb[sel]; A_list.append(np.full(len(b_sel), a, dtype=np.int64)); B_list.append(b_sel)
a = np.concatenate(A_list); b = np.concatenate(B_list); npts = len(a)
key = (a + B) * W + (b + B); order = np.argsort(key); key = key[order]; a = a[order]; b = b[order]
zabs = np.hypot(a, b)
# seed at 2+i: it is joined to 1+i by the unit step (1,0), so for every k >= 1 the component of 1+i is the
# component of 2+i plus the four ramified points; using 2+i lets us drop all odd-sum offsets (they can only
# touch 1+i and its associates, since every other Gaussian prime has a+b odd).
origin = int(np.searchsorted(key, (2 + B) * W + (1 + B))); assert a[origin] == 2 and b[origin] == 1
primes_count = int(isp.sum()); n_split = int(np.sum(np.nonzero(isp)[0] % 4 == 1)); n_inert_B = int(np.sum((np.nonzero(isp)[0] % 4 == 3) & (np.nonzero(isp)[0] <= B)))
assert npts == 8 * n_split + 4 * n_inert_B + 4, (npts, 8 * n_split + 4 * n_inert_B + 4)
R = {"radius": B, "norm_max": X, "gaussian_primes": npts, "sieve_s": round(t_sieve, 1), "points_s": round(time.time() - t0, 1)}
def component(K):
    rows = []; cols = []; r = math.isqrt(K)
    for dx in range(0, r + 1):
        for dy in range(-r, r + 1):
            if dx * dx + dy * dy > K or (dx == 0 and dy <= 0) or (dx + dy) % 2: continue   # odd dx+dy never joins two odd-parity primes
            tgt = key + dx * W + dy; pos = np.searchsorted(key, tgt); pos[pos >= npts] = 0
            hit = key[pos] == tgt; rows.append(np.nonzero(hit)[0].astype(np.int32)); cols.append(pos[hit].astype(np.int32)); del tgt, pos, hit
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    G = sparse.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(npts, npts))
    ncomp, lab = connected_components(G, directed=False); del G, rows, cols
    sel = lab == lab[origin]; reach = float(zabs[sel].max())
    return {"k": round(math.sqrt(K), 3), "k^2": K, "component_size": int(sel.sum()), "component_size_incl_1+i_associates": int(sel.sum()) + 4, "farthest_|z|": round(reach, 1),
            "farthest_norm": int(round(reach ** 2)), "closed_inside_disk": bool(reach <= B - math.sqrt(K)),
            "fraction_of_disk_primes": round(float(sel.mean()), 4),  "t_s": round(time.time() - t0, 1)}
chk = component(2); assert chk["component_size"] == 96, chk      # earlier receipt: 100 incl. the 4 ramified points
R["check_K2_component_of_2+i"] = chk["component_size"]
R["ladder"] = []
for K in (10, 16, 18, 20, 26, 32, 34, 36):
    info = component(K); R["ladder"].append(info); print(info, flush=True)
    if not info["closed_inside_disk"]: break
R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open(f"moat_rung2_receipts_R{B}.json", "w"), indent=1)
print({k: v for k, v in R.items() if k != "ladder"})
