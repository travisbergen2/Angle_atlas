"""angle_dial.py — The Angle Atlas on one page (2026-09-10).

Every angle is recomputed from its defining formula and asserted against the
thread records BEFORE drawing; nothing is typed in from memory.
Sources: Riemann-agent threads 2026-08-26 -> 2026-09-09 (see the Angle & Length
Atlas document). gamma_1 and gamma_100000 from Odlyzko's zeros1 table.
Nothing here bears on the truth of the Riemann Hypothesis.
Run: python3 angle_dial.py  -> angle_dial.png + angle_dial_receipts.txt
"""
import math, textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Arc, Circle

deg, rad = math.degrees, math.radians
G1 = 14.134725141734693        # gamma_1
G100K = 74920.827498994        # gamma_100000 (top of Odlyzko zeros1)

def rim(g):
    """Bergen Disc rim angle of a zero at height g under the Li map z = 1 - 1/s."""
    return deg(math.atan2(g, g * g - 0.25))

ang = {}
ang["rim_span"] = rim(G1) - rim(G100K)
ang["ladder"] = [90 / 2 ** k for k in range(1, 5)]            # 45, 22.5, 11.25, 5.625
ang["n45"] = deg(math.atan(11 / 45))                             # the number 45 on the sopfr dial
ang["atan_half"] = deg(math.atan(0.5))                           # arg(2+i); envelope slope 1/2
ang["bal35"] = deg(math.atan(0.35 / 0.65)); ang["bal35c"] = 90 - ang["bal35"]
ang["t345"] = deg(math.atan(4 / 3))                              # the 3-4-5 angle
ang["drift"] = [deg(math.atan((p + 1) / p)) for p in (2, 3, 5, 7, 11)]
ang["argrho1"] = deg(math.atan2(G1, 0.5)); ang["argrho100k"] = deg(math.atan2(G100K, 0.5))
ang["obtuse"] = deg(math.acos((4 + 9 - 16) / (2 * 2 * 3)))      # (2,3,4) largest angle
ang["phantom"] = deg(math.atan2(30, 0.2))                        # bearing of 0.7+30i from 1/2

# --- asserts against the thread records -------------------------------------
assert abs(ang["rim_span"] - 4.05) < 0.01
assert abs(ang["n45"] - 13.736) < 1e-3
assert abs(ang["atan_half"] - 26.565) < 1e-3
assert abs(ang["bal35"] - 28.30) < 0.01 and abs(ang["bal35c"] - 61.70) < 0.01
assert abs(ang["t345"] - deg(math.acos(3 / 5))) < 1e-9 and abs(ang["t345"] - 53.13) < 0.01
assert [round(x, 2) for x in ang["drift"][:3]] == [56.31, 53.13, 50.19]
assert abs(ang["argrho1"] - 87.97) < 0.01 and abs(ang["argrho100k"] - 89.9996) < 1e-4
assert abs(ang["obtuse"] - 104.48) < 0.01
assert abs(ang["phantom"] - 89.618) < 1e-3
assert 1 * 1 + 1 * (-1) == 0                                     # (1,1) . (1,-1): ramified pair is perpendicular
assert ang["ladder"] == [45.0, 22.5, 11.25, 5.625]

# --- colours -------------------------------------------------------------------
T, DEF, OPEN, KILL, MEAS, PAIR = "#1d4e89", "#2a9d8f", "#d62828", "#8d8d8d", "#7b2cbf", "#f77f00"

fig = plt.figure(figsize=(8.5, 11), dpi=200)
fig.patch.set_facecolor("white")
fig.text(0.5, 0.965, "The Angle Atlas — one page", ha="center", va="center", fontsize=16, fontweight="bold")
fig.text(0.5, 0.944, "Every load-bearing angle in the program on one dial · each value recomputed from its formula and asserted · 2026-09-10",
         ha="center", va="center", fontsize=8.5, color="#333")

axd = fig.add_axes([0.02, 0.50, 0.96, 0.42]); axd.set_aspect("equal"); axd.axis("off")
axd.set_xlim(-1.32, 1.32); axd.set_ylim(-0.12, 1.32)

def polar(a, r): return r * math.cos(rad(a)), r * math.sin(rad(a))
def ray(a, r, color, lw=2.2, ls="-", z=3):
    x, y = polar(a, r); axd.plot([0, x], [0, y], color=color, lw=lw, ls=ls, solid_capstyle="round", zorder=z)
def marker(a, r, n, color):
    x, y = polar(a, r)
    axd.add_patch(Circle((x, y), 0.044, facecolor=color, edgecolor="white", lw=1.2, zorder=6))
    axd.text(x, y, str(n), ha="center", va="center", fontsize=8.5, color="white", fontweight="bold", zorder=7)
def ring(a, r=1.0):
    x, y = polar(a, r); axd.add_patch(Circle((x, y), 0.058, facecolor="none", edgecolor=OPEN, lw=2.4, zorder=5))
def dlabel(a, r, text, color="#222", fs=7.5, bold=False):
    x, y = polar(a, r)
    axd.text(x, y, text, ha="center", va="center", fontsize=fs, color=color, zorder=8,
             fontweight="bold" if bold else "normal")

# rim, baseline, ticks
axd.add_patch(Arc((0, 0), 2, 2, theta1=0, theta2=180, lw=1.4, color="#444", zorder=2))
axd.plot([-1, 1], [0, 0], color="#444", lw=1.4, zorder=2)
for t in range(0, 181, 15):
    x0, y0 = polar(t, 1.0); x1, y1 = polar(t, 1.035 if t % 45 else 1.055)
    axd.plot([x0, x1], [y0, y1], color="#666", lw=0.9, zorder=2)
    if t in (30, 60, 120, 150):
        dlabel(t, 1.10, f"{t}°", color="#777", fs=7)
axd.add_patch(Circle((0, 0), 0.022, color="#222", zorder=9))

# 1  0°: real-axis mirror / inert primes / Weil-leg angle on the line
ray(0, 1.0, T, lw=2.6); ring(0); marker(0, 1.11, 1, T); dlabel(0, 1.225, "0°", fs=8, bold=True)
# 2  Disc rim span of the first 100,000 zeros
axd.add_patch(Wedge((0, 0), 0.5, 0, ang["rim_span"], facecolor=MEAS, alpha=0.45, edgecolor=MEAS, lw=0.8, zorder=2))
marker(2.0, 0.60, 2, MEAS)
# 3  the halving ladder 22.5 / 11.25 / 5.625
for a, r in zip(ang["ladder"][1:], (0.82, 0.72, 0.62)):
    ray(a, r, DEF, lw=2.0)
marker(22.5, 0.905, 3, DEF)
dlabel(11.25, 0.815, "11.25°", color=DEF, fs=6.5); dlabel(5.625, 0.72, "5.625°", color=DEF, fs=6.5)
# 4  the number 45 on the sopfr dial (killed '137')
ray(ang["n45"], 1.0, KILL, lw=1.6, ls="--"); marker(ang["n45"], 1.11, 4, KILL); dlabel(ang["n45"], 1.225, "13.7°", color=KILL, fs=7)
# 5  arctan 1/2
ray(ang["atan_half"], 1.0, T); marker(ang["atan_half"], 1.11, 5, T); dlabel(ang["atan_half"], 1.235, "26.57°", fs=7.5)
# 6  planted balance pair 28.30 / 61.70
ray(ang["bal35"], 0.62, PAIR, lw=1.8, ls="--"); ray(ang["bal35c"], 0.62, PAIR, lw=1.8, ls="--"); marker(ang["bal35c"], 0.70, 6, PAIR)
# 7  45°
ray(45, 1.0, T, lw=3.0); ring(45); marker(45, 1.11, 7, T); dlabel(45, 1.235, "45°", fs=8.5, bold=True)
# 8  the 3-4-5 angle
ray(ang["t345"], 0.90, T); marker(ang["t345"], 0.985, 8, T); dlabel(ang["t345"], 1.075, "53.13°", fs=7)
# 9  the (p, p+1) drift, ticks + bracket
for a in ang["drift"]:
    x0, y0 = polar(a, 0.965); x1, y1 = polar(a, 1.0); axd.plot([x0, x1], [y0, y1], color=T, lw=1.2, zorder=4)
axd.add_patch(Arc((0, 0), 2 * 1.165, 2 * 1.165, theta1=min(ang["drift"]), theta2=max(ang["drift"]), lw=1.2, color=T, zorder=4))
marker((min(ang["drift"]) + max(ang["drift"])) / 2, 1.245, 9, T)
# 10 arg rho wedge 87.97 -> 90
axd.add_patch(Wedge((0, 0), 0.80, ang["argrho1"], 90, facecolor=MEAS, alpha=0.45, edgecolor=MEAS, lw=0.8, zorder=2))
marker(84, 0.55, 10, MEAS)
# 11 90°
ray(90, 1.0, DEF, lw=3.0); ring(90); marker(90, 1.11, 11, DEF); dlabel(90, 1.235, "90°", fs=8.5, bold=True)
# 12 (2,3,4) obtuse
ray(ang["obtuse"], 0.62, KILL, lw=1.6, ls="--"); marker(ang["obtuse"], 0.70, 12, KILL); dlabel(ang["obtuse"], 0.82, "104.5°", color=KILL, fs=7)
# 13 E-RACE-A sector centers (0,45,90,135,180): tick at 135 + marker
x0, y0 = polar(135, 1.0); x1, y1 = polar(135, 1.06); axd.plot([x0, x1], [y0, y1], color=KILL, lw=1.6, zorder=4)
marker(135, 1.13, 13, KILL); dlabel(135, 1.25, "135°", color=KILL, fs=7)
# 14 180°: the pole factor's angle on the line
ray(180, 1.0, T, lw=2.4); marker(180, 1.11, 14, T); dlabel(180, 1.225, "180°", fs=8, bold=True)

# --- legend -------------------------------------------------------------------
axl = fig.add_axes([0.05, 0.02, 0.90, 0.46]); axl.axis("off"); axl.set_xlim(0, 1); axl.set_ylim(0, 1)
key = [(0.00, T, "theorem"), (0.13, DEF, "construction / definition"), (0.41, MEAS, "measured"),
       (0.55, PAIR, "planted exhibit"), (0.75, KILL, "killed")]
for x, col, name in key:
    axl.scatter([x + 0.008], [0.99], s=70, c=col, transform=axl.transAxes, clip_on=False)
    axl.text(x + 0.022, 0.99, name, fontsize=7.6, va="center", transform=axl.transAxes)
axl.scatter([0.008], [0.955], s=90, facecolors="none", edgecolors=OPEN, linewidths=2, transform=axl.transAxes, clip_on=False)
axl.text(0.022, 0.955, "red ring: the angle IS the Riemann Hypothesis at this value — at every zero (0°, 45°, 90°)", fontsize=7.6, va="center", color=OPEN, transform=axl.transAxes)

rows = [
 (1, "0°", "Real-axis mirror (conjugation); inert primes on the axis; the Weil leg's angle on the line.", "T · 0° ∀ zeros ⟺ RH", T),
 (2, f"{ang['rim_span']:.2f}°", "Bergen Disc rim occupied by all 100,000 known zeros, per side (rim angle ≈ 1/γ; #100000 sits 0.0053 px from the vanishing point).", "measured ✓", MEAS),
 (3, "22.5° · 11.25° · 5.625°", "The halving dial 90°/2^Ω(n) at Ω = 2, 3, 4 (27 sits at 11.25°); as roots of unity ζ₁₆, ζ₃₂, ζ₆₄ — one square root per rung.", "definition", DEF),
 (4, f"{ang['n45']:.2f}°", "The number 45 on the sopfr dial (arctan 11/45); the ‘137’ match (α⁻¹/10 = 13.704) dies at the third digit.", "KILLED", KILL),
 (5, f"{ang['atan_half']:.2f}°", "arctan ½ = arg(2+i); the composite envelope slope: sopfr(n) ≤ n/2 + 2, equality exactly on {4, 8} ∪ {2p}.", "T ✓", T),
 (6, f"{ang['bal35']:.1f}° / {ang['bal35c']:.1f}°", "Balance angles arctan(β/(1−β)) of a planted zero pair at β = 0.35: complementary (sum 90°, free) but not equal.", "T ✓ (exhibit)", PAIR),
 (7, "45°", "Primes on both dials (90°/2^Ω; arctan(sopfr(n)/n) = 45° ⟺ prime or 4); 1+i; a zero's balance angle is 45° ⟺ β = ½.", "T ✓ · 45° ∀ zeros ⟺ RH", T),
 (8, f"{ang['t345']:.2f}°", "The 3-4-5 angle = arg(3+4i) = arg((2+i)²); also the angle between 2+i and 2−i (arccos 3/5) — split factors are not perpendicular.", "T ✓", T),
 (9, f"{ang['drift'][0]:.1f}° → {ang['drift'][-1]:.1f}°", "arctan((p+1)/p) for p = 2, 3, 5, 7, 11: the (p, p+1) legs drift toward 45° and are never constant (the AAA argument's break 1).", "T ✓", T),
 (10, f"{ang['argrho1']:.2f}° → {ang['argrho100k']:.4f}°", "arg ρ — the zero waves' clock phase — from γ₁ to zero #100000.", "measured ✓", MEAS),
 (11, "90°", "Critical-line mirror s ↦ 1−s; the ramified pair 1±i is perpendicular; i is a quarter-turn; Paper 4's exponent 2; the duality vertex.", "T/defn · 90° ∀ zeros ⟺ RH", DEF),
 (12, f"{ang['obtuse']:.1f}°", "(2, 3, 4), the smallest scalene integer triangle, is obtuse — not a right triangle.", "KILLED (as ‘right’)", KILL),
 (13, "0°, 45°, 90°, 135°, 180°", "E-RACE-A sector centers (half-widths 22.5° to 180°); the only wall-free sharp sector is S(0°, ±22.5°), the axis-hugging window.", "measured ✓", KILL),
 (14, "180°", "The pole factor's angle on the line; pole + Gamma + ζ angles ≡ 0 mod 180° (Z(t) is real).", "T ✓", T),
]
y, dy = 0.905, 0.060
for n, a, what, grade, col in rows:
    axl.scatter([0.012], [y], s=150, c=col, transform=axl.transAxes, clip_on=False, zorder=3)
    axl.text(0.012, y, str(n), fontsize=7, color="white", fontweight="bold", ha="center", va="center", transform=axl.transAxes, zorder=4)
    axl.text(0.036, y, a, fontsize=7.8, fontweight="bold", va="center", transform=axl.transAxes)
    lines = textwrap.wrap(what, 72)
    axl.text(0.245, y, "\n".join(lines), fontsize=7.4, va="center", linespacing=1.12, transform=axl.transAxes)
    axl.text(0.815, y, grade, fontsize=7.0, color=col, va="center", transform=axl.transAxes)
    y -= dy

fig.text(0.5, 0.008, "Red rings mark the three places where the angle IS the Riemann Hypothesis (0° · 45° · 90°) — at every zero; verified for 100,000, proven for none. "
         "Nothing here bears on RH's truth.\nSources: Riemann-agent threads 2026-08-26 → 09-09 · γ₁, γ₁₀₀₀₀₀ from Odlyzko · script angle_dial.py (asserts pass)",
         ha="center", va="bottom", fontsize=6.6, color="#555", linespacing=1.3)

out = "/agent/workspace/angle_atlas/angle_dial.png"
fig.savefig(out, dpi=200, facecolor="white")
with open("/agent/workspace/angle_atlas/angle_dial_receipts.txt", "w") as f:
    for k, v in ang.items(): f.write(f"{k}: {v}\n")
    f.write("all asserts passed\n")
print("saved", out); print({k: (round(v, 4) if isinstance(v, float) else [round(x, 4) for x in v]) for k, v in ang.items()})
