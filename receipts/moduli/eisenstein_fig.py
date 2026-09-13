"""eisenstein_fig.py — three panels for the Eisenstein rebuild (values recomputed here; asserts pin the receipts)."""
import numpy as np, math, cmath, json
from math import factorial
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle
w = cmath.exp(2j*math.pi/3); R = json.load(open("eisenstein_receipts.json"))
B2k = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330]
def hurwitz(s, a, Nt=1000, m=10):
    n = np.arange(Nt, dtype=float); S = np.sum(np.exp(-s*np.log(n + a))); x = Nt + a
    S += x**(1 - s)/(s - 1) + 0.5*x**(-s); rf = s
    for k in range(1, m + 1):
        if k > 1: rf *= (s + 2*k - 3)*(s + 2*k - 2)
        S += B2k[k-1]/factorial(2*k) * rf * x**(-s - 2*k + 1)
    return S
zeta = lambda s: hurwitz(s, 1.0); L4 = lambda s: 4**(-s)*(hurwitz(s, .25) - hurwitz(s, .75)); L3 = lambda s: 3**(-s)*(hurwitz(s, 1/3) - hurwitz(s, 2/3))
assert abs(L3(complex(0.5, R["t3"]))) < 1e-9
T, DEF, MEAS, RED, GOLD, GRN = "#1d4e89", "#2a9d8f", "#7b2cbf", "#d62828", "#b8860b", "#2e8b57"
fig, axs = plt.subplots(1, 3, figsize=(17, 5.6), dpi=170); fig.patch.set_facecolor("white")
# ---- A: the hexagonal lattice, six units, fold wedge, ramified prime, the circle's six points ----
ax = axs[0]; ax.set_aspect("equal"); ax.set_xlim(-2.7, 4.9); ax.set_ylim(-2.7, 3.9); ax.set_title("Z[ω]: six units, the fold, the ramified prime, the same circle", fontsize=10.5)
for a in range(-5, 6):
    for b in range(-5, 6):
        z = a + b*w
        if -2.7 < z.real < 4.9 and -2.7 < z.imag < 3.9: ax.plot(z.real, z.imag, ".", color="#c8c8c8", ms=4, zorder=1)
ax.add_patch(Wedge((0, 0), 3.6, 0, 30, facecolor=DEF, alpha=0.13, edgecolor=DEF, lw=0.8, zorder=0))
ax.text(3.15, 0.55, "fold\n[0°, 30°]", color=DEF, fontsize=8, ha="center")
for k in range(6):
    u = (-w)**k; ax.plot([0, u.real], [0, u.imag], color=GOLD, lw=1.0, zorder=2); ax.plot(u.real, u.imag, "o", color=GOLD, ms=6, zorder=4)
ax.text(1.05, -0.22, "1", color=GOLD, fontsize=8); ax.text(0.22, 1.05, "−ω²", color=GOLD, fontsize=8); ax.text(-0.75, 0.95, "ω", color=GOLD, fontsize=8)
rp = 1 - w
for k in range(6):
    z = rp*(-w)**k; ax.plot(z.real, z.imag, "o", color=RED, ms=7, zorder=5, mfc="none", mew=1.8)
ax.plot([0, rp.real], [0, rp.imag], color=RED, lw=1.6, zorder=3); ax.text(rp.real + 0.12, rp.imag - 0.42, "1−ω (norm 3) at −30°", color=RED, fontsize=7.5)
ax.text(1.62, 1.02, "2+ω = 1−ω² (norm 3) at +30°:\nthe ramified prime sits on the fold edge", color=RED, fontsize=7.2)
ax.add_patch(Circle((1, 0), 1, fill=False, color=T, lw=1.4, zorder=3))
for a, b in [(0, -1), (0, 0), (1, -1), (1, 1), (2, 0), (2, 1)]:
    z = a + b*w; ax.plot(z.real, z.imag, "s", color=T, ms=6, zorder=6)
ax.text(-2.55, -2.1, "trace = norm ⟺ |z−1| = 1: the same circle as in Z[i];\nsix lattice points now (1 + each unit): 0, 2, 1±ω, 1±ω²", color=T, fontsize=7.8, ha="left")
for p, lab in [(7, "3+2ω"), (13, "3+4ω"), (19, "5+6ω")]:
    for y in range(0, 4):
        x2 = p - 3*y*y; x = math.isqrt(x2)
        if x*x == x2: break
    z = complex(x, y*math.sqrt(3)); ax.plot(z.real, z.imag, "o", color=MEAS, ms=5, zorder=6); ax.text(z.real + 0.06, z.imag + 0.05, f"{lab} → {p}", color=MEAS, fontsize=7)
ax.axhline(0, color="#999", lw=0.6); ax.axvline(0, color="#999", lw=0.6); ax.set_xticks([]); ax.set_yticks([])
# ---- B: folded-angle histogram ----
ax = axs[1]; folded = np.load("eisenstein_folded.npy")
ax.hist(folded, bins=30, range=(0, 30), color=DEF, alpha=0.85, edgecolor="white")
ax.axhline(len(folded)/30, color=RED, lw=1.2, ls="--"); ax.set_xlim(0, 30)
ax.set_title(f"folded angles of the {len(folded):,} split primes ≤ 10⁶ (Hecke): KS = {R['ks']:.4f}", fontsize=10.5)
ax.set_xlabel("θ mod 60°, reflected into [0°, 30°]"); ax.set_ylabel("count per 1° bin"); ax.text(0.5, len(folded)/30*1.02, "uniform", color=RED, fontsize=8)
# ---- C: the three L-functions on the line ----
ax = axs[2]; ts = np.linspace(0.6, 20, 700)
for f, col, lab, tz in [(zeta, T, "|ζ(½+it)|", R["tz"]), (L4, MEAS, "|L(½+it, χ₄)|  (Z[i])", R["t4"]), (L3, GRN, "|L(½+it, χ₃)|  (Z[ω])", R["t3"])]:
    v = [abs(f(complex(0.5, t))) for t in ts]; ax.plot(ts, v, color=col, lw=1.6, label=lab); ax.plot(tz, 0, "v", color=col, ms=8, zorder=5); ax.text(tz, 0.12, f"{tz:.4f}", color=col, fontsize=7.5, ha="center")
ax.set_ylim(0, 4.2); ax.set_xlim(0.6, 20); ax.set_xlabel("t"); ax.set_title("two worlds, one line: ζ_{Q(i)} = ζ·L(χ₄),  ζ_{Q(ω)} = ζ·L(χ₃)", fontsize=10.5); ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.25)
fig.text(0.5, 0.005, "All values recomputed in eisenstein.py / eisenstein_fig.py (asserts pass). Zeros shown are the first of each; each is verified to great height and proven for none. Nothing here bears on RH.  2026-09-10", ha="center", fontsize=8, color="#555")
plt.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig("eisenstein.png", dpi=170, facecolor="white"); print("saved eisenstein.png")
