"""reciprocal_check.py — receipts for 'the reciprocals of the zeta zeros are the other end of the line' (2026-09-10).
(1) Under w = 1/s the critical line Re s = 1/2 maps to the circle |w - 1| = 1 (through 0 = image of infinity, and 2 = image of s = 1/2).
(2) That circle is the trace = norm circle (a-1)^2 + b^2 = 1 of the 09-04 Warren Note: trace(w) = norm(w) <=> 2 Re w = |w|^2 <=> Re(1/w) = 1/2.
(3) On the Bergen Disc (Li map z = 1 - 1/s): 1 - z_rho = 1/rho, so 1/rho is the vector from the vanishing point to the zero's image; |z_rho| = 1 <=> 1/rho on the circle <=> beta = 1/2.
(4) The reciprocals of all zeros sit within 1/|rho_1| = 0.0707 of the origin — the far end of the line, seen up close.
"""
import math, cmath
# Points on the critical line. The identity needs only beta = 1/2; heights are the first three zeros
# (Odlyzko table values, as used in verify_angles.py) plus an arbitrary height, to show no zero data is needed.
HEIGHTS = [14.134725141734693, 21.022039638771555, 25.010857580145688, 1000.5]
zeros = [complex(0.5, g) for g in HEIGHTS]
src = "beta = 1/2 at heights " + ", ".join(f"{g:g}" for g in HEIGHTS) + " (identity needs only beta = 1/2)"
print("zero source:", src)
if zeros:
    g1 = zeros[0].imag
    print(f"gamma_1 = {g1:.12f}")
    worst = max(abs(abs(1/r - 1) - 1) for r in zeros)
    print(f"(1) max | |1/rho - 1| - 1 | over the on-line points = {worst:.2e}  (0 => reciprocals lie on the circle |w-1|=1)")
    assert worst < 1e-12
    print(f"(4) |1/rho_1| = {abs(1/zeros[0]):.6f}; |1/(1/2+1000.5i)| = {abs(1/zeros[3]):.6f}; all reciprocals of on-line points with height >= gamma_1 lie within {max(abs(1/r) for r in zeros):.4f} of the origin")
    z1 = 1 - 1/zeros[0]
    print(f"(3) Li map: z_rho1 = {z1:.6f}, |z_rho1| = {abs(z1):.12f}; 1 - z = 1/rho: {abs((1 - z1) - 1/zeros[0]):.1e}")
# planted off-line zero
rho_off = complex(0.6, 30.0)
w = 1/rho_off
print(f"planted beta=0.6: |1/rho - 1| = {abs(w-1):.6f} (not 1); |z| = |1 - 1/rho| = {abs(1 - w):.6f}; predicted |z|^2 - 1 = (1-2beta)/|rho|^2 = {(1-1.2)/abs(rho_off)**2:.6e}, actual {abs(1-w)**2-1:.6e}")
# (2) the identity, on the four Gaussian lattice points of the trace = norm circle and on random points
for name, wv in [("0", 0), ("2", 2), ("1+i", complex(1, 1)), ("1-i", complex(1, -1))]:
    if wv == 0:
        print(f"(2) w = 0: trace 0 = norm 0; reciprocal = the point at infinity of the line")
    else:
        tr, nm = 2*wv.real if isinstance(wv, complex) else 2*wv, abs(wv)**2
        rec = 1/wv
        print(f"(2) w = {name}: trace {tr:g} = norm {nm:g}; 1/w = {rec}  Re = {rec.real:g}  -> on the critical line at height {(-rec.imag if isinstance(rec, complex) else 0):g}")
import random
random.seed(1)
ok = 0
for _ in range(10000):
    wv = complex(random.uniform(-3, 3), random.uniform(-3, 3))
    if abs(wv) < 1e-9: continue
    lhs = abs(2*wv.real - abs(wv)**2) < 1e-9
    rhs = abs((1/wv).real - 0.5) < 1e-9
    # test the identity as an equation: 2Re w - |w|^2 = |w|^2 (2 Re(1/w) - 1) exactly
    ok += abs((2*wv.real - abs(wv)**2) - abs(wv)**2*(2*(1/wv).real - 1)) < 1e-9
print(f"(2) identity 2Re w - |w|^2 = |w|^2 (2Re(1/w) - 1) holds on {ok}/10000 random points")
# which Gaussian integers have reciprocal on the line? enumerate |a|,|b| <= 60
pts = [(a, b) for a in range(-60, 61) for b in range(-60, 61) if (a, b) != (0, 0) and abs(complex(a, b) and (1/complex(a, b)).real - 0.5) < 1e-12]
print("(2) Gaussian integers (|a|,|b|<=60) whose reciprocal lies on Re = 1/2:", pts)
assert sorted(pts) == [(1, -1), (1, 1), (2, 0)]
# heights on the line of those reciprocals, vs the first zero
print("    their reciprocals sit at heights 0 and ±1/2 on the line; the first zero is at height 14.13 — the zeros' reciprocals crowd the far end (near 0), the lattice points sit at the near end (near 2).")
print("ALL CHECKS PASSED")
