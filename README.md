# The Angle Atlas — the dial

**Live:** https://fractalyouniverse.org/Angle_atlas/ · Room VII of [Fractal Youniverse](https://fractalyouniverse.org/)

Every load-bearing angle in the Information Manifold Model's Riemann program, on one half-dial (0° to 180°). Tap a ray, or a numbered marker, or a legend row, to read what the angle is, where it comes from, and its grade. Arrow keys and the ‹ › buttons step through the fourteen rays.

**Instruments, not proofs.** Nothing on this page bears on the truth of the Riemann Hypothesis.

## What is on the dial

| # | angle | what | grade |
|---|---|---|---|
| 1 | 0° | the real-axis mirror; inert primes; the Weil leg's angle on the line | theorem · **RH lives here** |
| 2 | 4.05° | the rim of the Bergen Disc occupied by all 100,000 tabulated zeros (per side) | measured |
| 3 | 22.5° · 11.25° · 5.625° | the halving dial 90°/2^Ω(n) at Ω = 2, 3, 4 | definition |
| 4 | 13.74° | the number 45 on the sopfr dial — the "137" coincidence, killed | killed |
| 5 | 26.57° | arctan ½ = arg(2+i); the composite envelope slope | theorem |
| 6 | 28.3° / 61.7° | balance angles of a planted zero pair at β = 0.35 — complementary, not equal | planted exhibit |
| 7 | 45° | primes on both dials; 1+i; a zero's balance angle iff β = ½ | theorem · **RH lives here** |
| 8 | 53.13° | the 3-4-5 angle = arg((2+i)²); the angle between 2+i and 2−i | theorem |
| 9 | 56.3° → 47.5° | arctan((p+1)/p): the (p, p+1) legs drift toward 45°, never constant | theorem |
| 10 | 87.97° → 89.9996° | arg ρ — the zero waves' clock phases | measured |
| 11 | 90° | the critical-line mirror; the ramified pair 1±i perpendicular; the bearing from ½ | theorem / construction · **RH lives here** |
| 12 | 104.5° | (2, 3, 4) is obtuse, not right | killed |
| 13 | 0°, 45°, 90°, 135°, 180° | the angular race census's sector centers; only S(0°, ±22.5°) is wall-free | measured |
| 14 | 180° | the pole factor's angle on the line; pole + Gamma + ζ ≡ 0 mod 180° | theorem |

**The three red rings** (rays 1, 7, 11) mark the places where the angle *is* the Riemann Hypothesis — Weil-leg angle 0°, balance angle 45°, bearing 90°, each at every zero. Each is β = ½ rewritten; each is verified for every zero ever computed and proven for none.

**The fence.** Every ray is one of three kinds — free (true for any function with the mirror and the strip), equivalent (RH rewritten), or measured (a finite window). The Davenport–Heilbronn function (1936) passes every free ray and has zeros off the line, so any argument built from this dial alone proves a falsehood; the missing ingredient is the Euler product over all primes.

## Receipts

Every angle is recomputed **in the browser** from its defining formula and compared with the values asserted by the committed Python twin (`receipts/angle_dial.py`); the page reports `15 / 15 recomputed values agree`. The Python side:

```bash
python3 receipts/verify_angles.py      # recomputes and asserts every quoted angle and length (no dependencies)
python3 receipts/angle_dial.py         # asserts, then draws the one-page figure angle_dial.png (needs matplotlib)
python3 receipts/reciprocal_check.py   # the reciprocal-of-the-line identity: Re(1/w) = ½ ⟺ trace(w) = norm(w)
```

Constants: γ₁ = 14.134725141734693 and γ₁₀₀₀₀₀ = 74920.827498994 from Odlyzko's `zeros1` table.

## Files

- `index.html` — the instrument, self-contained (one external request: the font stylesheet; no analytics, no tracking)
- `receipts/angle_dial.py` — Python twin: asserts every angle, draws `angle_dial.png`
- `receipts/verify_angles.py` — the atlas's recomputation script
- `receipts/reciprocal_check.py` — the reciprocal-circle identity
- `angle_dial.png` — the one-page figure

## Embed

```html
<iframe src="https://fractalyouniverse.org/Angle_atlas/" width="100%" height="900" style="border:0;border-radius:14px" loading="lazy" title="The Angle Atlas"></iframe>
```

## Sources

- T. Bergen, *A Smith Chart for the Riemann Hypothesis* (IMM Paper 19), doi:10.5281/zenodo.21400881
- The companion instruments: [The RH Smith Chart](https://travisbergen2.github.io/rh-smith-chart/) and [The Bergen Disc](https://fractalyouniverse.org/Bergen_disc/)
- A. M. Odlyzko, tables of zeros of the Riemann zeta function
- E. Hecke, Math. Z. 1 (1918) 357–376; 6 (1920) 11–51
- H. Davenport, H. Heilbronn, J. London Math. Soc. s1-11 (1936) 181–185, 307–312
- D. Platt, T. Trudgian, Bull. London Math. Soc. 53 (2021) 792–797

Built 2026-09-10 by Travis Bergen with the Riemann agent.
