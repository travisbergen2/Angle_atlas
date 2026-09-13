"""gaussian_moat.py — the Gaussian moat visual (Travis, 2026-09-11 09:27).

Gaussian primes with norm <= 1e6 (disk of radius 1000), split (norm p) in one
colour, inert (on the axes, norm p^2) in another; the connected component of
1+i under steps of length <= k for a ladder of k, each bounded by its moat;
and the largest prime-free circle inside the disk — the 2-D record gap.
Everything computed here; the figure is drawn only after the asserts pass.
"""
import numpy as np, math, json, time
from scipy import ndimage, sparse
from scipy.sparse.csgraph import connected_components
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

t0 = time.time(); R = {}
B = 2000; X = B * B          # disk of radius 2000, norms to 4e6
# ---- rational primes to 1e6 ------------------------------------------------
isp = np.ones(X + 1, dtype=bool); isp[:2] = False
for p in range(2, int(X**0.5) + 1):
    if isp[p]: isp[p*p::p] = False
primes = np.nonzero(isp)[0]
n_split = int(np.sum(primes % 4 == 1)); n_inert_le_B = int(np.sum((primes % 4 == 3) & (primes <= B)))

# ---- Gaussian primes on the lattice -----------------------------------------
a = np.arange(-B, B + 1, dtype=np.int32); A, Bm = np.meshgrid(a, a, indexing="ij")     # A = real part, Bm = imag part
nm = A * A + Bm * Bm
split = (nm <= X) & isp[np.minimum(nm, X)]                             # a^2+b^2 prime  (=> 2 or 1 mod 4)
inert = ((A == 0) & isp[np.abs(Bm)] & (np.abs(Bm) % 4 == 3)) | ((Bm == 0) & isp[np.abs(A)] & (np.abs(A) % 4 == 3))
M = split | inert
R["gaussian_primes_norm<=1e6"] = int(M.sum())
R["expected_8*split+4*inert(p<=1000)+4"] = 8 * n_split + 4 * n_inert_le_B + 4
assert R["gaussian_primes_norm<=1e6"] == R["expected_8*split+4*inert(p<=1000)+4"]
R["split_points"] = int(split.sum()) - 4; R["inert_points"] = int(inert.sum()); R["ramified_points(1+i and associates)"] = 4
assert bool(M[B + 1, B + 1])                                             # 1+i is there
dens_center = float(M[B-50:B+51, B-50:B+51].mean()); dens_rim = float(M[(nm > 0.9*X) & (nm <= X)].mean())
R["density_center_101x101"] = round(dens_center, 4); R["density_rim_shell_0.9-1.0"] = round(dens_rim, 4)
R["(4/pi)/log_N_at_rim"] = round((4 / math.pi) / math.log(0.95 * X), 4)

# ---- components of 1+i under steps <= k ----------------------------------------
idx = -np.ones(M.shape, dtype=np.int32); pts = np.argwhere(M); idx[M] = np.arange(len(pts))
npts = len(pts); zabs = np.hypot(pts[:, 0] - B, pts[:, 1] - B)
def component(K):
    rows = []; cols = []
    for dx in range(0, int(math.isqrt(K)) + 1):
        for dy in range(-int(math.isqrt(K)), int(math.isqrt(K)) + 1):
            if dx * dx + dy * dy > K or (dx == 0 and dy <= 0): continue
            s1 = (slice(max(0, -dx), M.shape[0] - max(0, dx)), slice(max(0, -dy), M.shape[1] - max(0, dy)))
            s2 = (slice(max(0, dx), M.shape[0] - max(0, -dx)), slice(max(0, dy), M.shape[1] - max(0, -dy)))
            both = M[s1] & M[s2]
            rows.append(idx[s1][both]); cols.append(idx[s2][both])
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    G = sparse.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(npts, npts))
    ncomp, lab = connected_components(G, directed=False)
    c = lab[idx[B + 1, B + 1]]; sel = lab == c
    reach = float(zabs[sel].max())
    return sel, {"k": round(math.sqrt(K), 3), "k^2": K, "component_size": int(sel.sum()), "farthest_|z|": round(reach, 1),
                 "farthest_norm": int(round(reach ** 2)), "closed_inside_disk": bool(reach <= B - math.sqrt(K)), "components_total": int(ncomp)}
# only EVEN k^2 can join two Gaussian primes other than 1+i: every other Gaussian prime has a+b odd,
# so a step (dx,dy) between two of them needs dx+dy even, i.e. dx^2+dy^2 even.  (Checked below.)
ladder = [2, 4, 8, 10, 16, 18, 20, 26, 32, 34, 36]
R['odd_k^2_useless_check'] = component(5)[1]['component_size'] == component(4)[1]['component_size']
comps = {}; rows_out = []
for K in ladder:
    sel, info = component(K); comps[K] = sel; rows_out.append(info)
    if not info["closed_inside_disk"]: break
R["moat_ladder"] = rows_out
closed = [r for r in rows_out if r["closed_inside_disk"]]
R["largest_k_closed_inside_radius_1000"] = closed[-1]["k"] if closed else None
R["first_k_reaching_boundary(undetermined_here)"] = rows_out[-1]["k"] if not rows_out[-1]["closed_inside_disk"] else None

# ---- the largest prime-free circle inside the disk (the 2-D record gap) ---------
D = ndimage.distance_transform_edt(~M)                 # distance from each non-prime point to nearest Gaussian prime
rad = np.hypot(A, Bm); fits = D <= (B - rad)           # circle must lie inside the disk
Dm = np.where(fits, D, 0.0); k_ = np.unravel_index(int(Dm.argmax()), Dm.shape)
r_max = float(Dm[k_]); ca, cb = int(A[k_]), int(Bm[k_])
inside = int(np.sum((A - ca) ** 2 + (Bm - cb) ** 2 < r_max ** 2))
R["largest_empty_circle"] = {"radius": round(r_max, 3), "center": (ca, cb), "center_norm": ca * ca + cb * cb,
                             "lattice_points_inside": inside, "primes_inside": int(np.sum(M & ((A - ca) ** 2 + (Bm - cb) ** 2 < r_max ** 2)))}
# Poisson-scale expectation for the largest empty circle near the rim: pi r^2 * density ~ log(#points)
R["poisson_scale_radius_estimate"] = round(math.sqrt(math.log(math.pi * X) / (math.pi * dens_rim)), 2)
pg = np.diff(primes); kk = int(pg.argmax())
R["compare_1D_record_gap_below_X"] = {"X": X, "gap": int(pg[kk]), "after_prime": int(primes[kk]), "composites_inside": int(pg[kk]) - 1}
R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open("gaussian_moat_receipts.json", "w"), indent=1, default=str)
for k, v in R.items():
    if k == "moat_ladder": print(k); [print("   ", r) for r in v]
    else: print(k, v)
assert R["largest_empty_circle"]["primes_inside"] == 0
print("ALL ASSERTS PASS — drawing")

# ---- figure ------------------------------------------------------------------------
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
fig = plt.figure(figsize=(16, 6.2), dpi=150)
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.22, top=0.84, bottom=0.08)
C_SPLIT = "#2b3f6b"; C_INERT = "#c8102e"; C_RAM = "#e8a300"
cmap = plt.cm.viridis
# A: the disk, components overlaid
ax = fig.add_subplot(gs[0, 0]); ax.set_aspect("equal"); ax.set_facecolor("white")
img = np.full(M.shape + (3,), 255, dtype=np.uint8); img[split] = (np.array(matplotlib.colors.to_rgb("#9fb0c9")) * 255).astype(np.uint8); img[inert] = (np.array(matplotlib.colors.to_rgb(C_INERT)) * 255).astype(np.uint8)
ax.imshow(np.transpose(img, (1, 0, 2)), origin="lower", extent=[-B - .5, B + .5, -B - .5, B + .5], interpolation="nearest")
ncl = len(closed)
COL = ["#3b0f70", "#1f77b4", "#1b9e77", "#f2c14e", "#e377c2", "#8c564b"]          # smallest k first
for i, r in reversed(list(enumerate(closed))):                                    # draw largest first so the small ones stay visible
    sel = comps[r["k^2"]]; P = pts[sel]
    ax.scatter(P[:, 0] - B, P[:, 1] - B, s=0.6 + 1.2 * (ncl - 1 - i), color=COL[i], lw=0,
               label=f"k = √{r['k^2']}: {r['component_size']:,} primes, moat beyond |z| = {r['farthest_|z|']:.0f}")
ax.add_patch(Circle((0, 0), closed[-1]["farthest_|z|"], fill=False, lw=0.9, ls="--", color="#333"))
ax.annotate(f"√{closed[-1]['k^2']} moat: no Gaussian prime of the component\nlies beyond |z| = {closed[-1]['farthest_|z|']:.0f}",
            xy=(-closed[-1]["farthest_|z|"] * 0.71, closed[-1]["farthest_|z|"] * 0.71), xytext=(-1950, 1560), fontsize=7,
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#333"))
ax.add_patch(Circle((0, 0), B, fill=False, lw=0.8, color="#333"))
ax.set_xlim(-B - 20, B + 20); ax.set_ylim(-B - 20, B + 20)
ax.set_title(f"Gaussian primes, norm ≤ {X:,}  ({R['gaussian_primes_norm<=1e6']:,} points)\n"
             f"grey: split (norm p)  ·  red: inert on the axes (norm p²)  ·  colour: the component of 1+i under steps ≤ k", fontsize=9)
ax.legend(loc="lower left", fontsize=6.5, markerscale=6, framealpha=0.95, title="components of 1+i closed inside the disk", title_fontsize=7)
ax.set_xlabel("Re"); ax.set_ylabel("Im")
# B: zoom near the origin with the nested components and their moats
ax = fig.add_subplot(gs[0, 1]); ax.set_aspect("equal"); Z = 60
sub = (slice(B - Z, B + Z + 1), slice(B - Z, B + Z + 1))
sp = np.argwhere(split[sub]); ip = np.argwhere(inert[sub])
ax.scatter(sp[:, 0] - Z, sp[:, 1] - Z, s=6, color="#9fb0c9", lw=0)
ax.scatter(ip[:, 0] - Z, ip[:, 1] - Z, s=8, color=C_INERT, lw=0)
for i, r in reversed(list(enumerate(closed[:6]))):
    P = pts[comps[r["k^2"]]] - B; P = P[(np.abs(P[:, 0]) <= Z) & (np.abs(P[:, 1]) <= Z)]
    ax.scatter(P[:, 0], P[:, 1], s=10 + 7 * (ncl - 1 - i), color=COL[i], lw=0, label=f"k = √{r['k^2']}  (reach {r['farthest_|z|']:.0f})")
ax.scatter([1, -1, 1, -1], [1, 1, -1, -1], s=30, color=C_RAM, marker="D", lw=0, label="1+i and associates")
ax.set_xlim(-Z, Z); ax.set_ylim(-Z, Z); ax.set_title("near the origin: the nested components of 1+i\n(larger markers = smaller step; each is stopped by its own moat)", fontsize=8.5)
ax.legend(loc="upper right", fontsize=6.5, framealpha=0.95); ax.set_xlabel("Re"); ax.set_ylabel("Im")
# C: the largest prime-free circle
ax = fig.add_subplot(gs[0, 2]); ax.set_aspect("equal"); W = int(r_max) + 14
sub = (slice(B + ca - W, B + ca + W + 1), slice(B + cb - W, B + cb + W + 1))
sp = np.argwhere(split[sub]); ip = np.argwhere(inert[sub])
ax.scatter(sp[:, 0] - W + ca, sp[:, 1] - W + cb, s=10, color=C_SPLIT, lw=0, label="Gaussian prime (split)")
if len(ip): ax.scatter(ip[:, 0] - W + ca, ip[:, 1] - W + cb, s=12, color=C_INERT, lw=0, label="inert")
ax.add_patch(Circle((ca, cb), r_max, fill=True, alpha=0.12, color="#c8102e", lw=0))
ax.add_patch(Circle((ca, cb), r_max, fill=False, lw=1.4, color="#c8102e"))
ax.plot([ca], [cb], "+", color="#c8102e", ms=9)
ax.set_xlim(ca - W, ca + W); ax.set_ylim(cb - W, cb + W)
ax.set_title(f"the largest prime-free circle in the disk: radius {r_max:.2f} at {ca}{cb:+d}i, norm {ca*ca+cb*cb:,}\n"
             f"{inside} lattice points, 0 primes — the 2-D record gap\n(1-D record below {X:,}: gap {R['compare_1D_record_gap_below_X']['gap']} after {R['compare_1D_record_gap_below_X']['after_prime']:,}, {R['compare_1D_record_gap_below_X']['composites_inside']} composites)", fontsize=8)
ax.legend(loc="upper right", fontsize=6.5, framealpha=0.95); ax.set_xlabel("Re"); ax.set_ylabel("Im")
fig.suptitle("The Gaussian moat — the other end of the line, in two dimensions.   Instruments, not proofs: whether moats of every width exist is open.", fontsize=10, y=0.985)
fig.savefig("gaussian_moat.png", bbox_inches="tight", facecolor="white")
print("figure written", round(time.time() - t0, 1), "s")
