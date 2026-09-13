"""trace_norm_rh.py — RH as 'trace = norm' for the reciprocal zeros (2026-09-10).
w = 1/rho. trace(w) = 2 Re w, norm(w) = |w|^2. Identity: trace(w) - norm(w) = |w|^2 (2 Re(1/w) - 1) = |w|^2 (2 beta - 1), beta = Re rho.
Mirror pair rho, 1 - conj(rho) (betas beta, 1-beta): sum of [norm - trace] = (2 beta - 1) (1/|1-conj rho|^2 - 1/|rho|^2) >= 0, = 0 iff beta = 1/2.
So D := sum over mirror pairs [norm(1/rho) - trace(1/rho)] >= 0 always, and D = 0 <=> RH.
Trace side: sum trace(1/rho) = 2 sum Re(1/rho) = 2 lambda_1 = 2 + gamma_E - log(4 pi) (unconditional, arithmetic).
Li's criterion (1997): lambda_n = sum_rho [1 - (1 - 1/rho)^n] >= 0 for all n  <=>  RH: positivity of trace data forces |1 - 1/rho| = 1."""
import math, cmath
gE = 0.5772156649015329
lam1 = 1 + gE/2 - 0.5*math.log(4*math.pi)
print(f"lambda_1 = sum_rho 1/rho = 1 + gamma/2 - (1/2) log 4pi = {lam1:.7f};  sum trace(1/rho) = 2 lambda_1 = {2*lam1:.7f}  [T, unconditional]")
print("mirror pair (beta, gamma): norm - trace summed over the pair, and the closed form (2b-1)(1/|1-conj rho|^2 - 1/|rho|^2):")
for beta, gamma in [(0.5, 14.134725), (0.6, 30.0), (0.7, 30.0), (0.9, 49.5), (0.5, 1000.5)]:
    rho = complex(beta, gamma); mir = 1 - rho.conjugate()
    def nt(r): w = 1/r; return abs(w)**2 - 2*w.real
    D = nt(rho) + nt(mir); closed = (2*beta - 1)*(1/abs(mir)**2 - 1/abs(rho)**2)
    assert abs(D - closed) < 1e-15 and D >= -1e-18
    print(f"   beta={beta:.1f}, gamma={gamma:>8.4f}: D_pair = {D:+.3e}  closed form {closed:+.3e}   {'= 0 (on the line)' if abs(D) < 1e-15 else '> 0 (off the line)'}")
print("   -> D >= 0 always; D = 0 iff every zero has beta = 1/2: RH is 'trace = norm' for the reciprocal zeros  [T]")
print("\nLi's mechanism on a toy pair z, 1/conj(z): Re sum (1 - z^n) over the pair = 2 - (|z|^n + |z|^-n) cos(n theta)")
for r in (1.0, 1.02):
    z = r*cmath.exp(1j*math.radians(20)); zm = 1/z.conjugate()
    vals = [(1 - z**n).real + (1 - zm**n).real for n in range(1, 400)]
    mn = min(vals); nmin = vals.index(mn) + 1
    print(f"   |z| = {r}: min over n<400 of the pair's Li contribution = {mn:+.4f} at n = {nmin}  -> {'stays >= 0' if mn >= -1e-12 else 'goes NEGATIVE: off the circle is detected by some n'}")
print("ALL CHECKS PASSED")
