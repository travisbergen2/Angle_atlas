"""product_edge_receipts.py — Python twin for Room VIII (The Product's Edge), 2026-09-10.
Computes every checkpoint the page recomputes in the browser and writes product_edge_receipts.json.
No mpmath here: zero heights are Odlyzko-table values, each VERIFIED below by |zeta(1/2 + i gamma_k)| ~ 0 via the
pair-averaged eta series (independent of any zero table)."""
import numpy as np, math, json
N = 10**6
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N**0.5) + 1):
    if spf[p] == 0:
        block = spf[p*p::p]; block[block == 0] = p
spf[spf == 0] = np.arange(N + 1)[spf == 0]; spf[0] = spf[1] = 0
n = np.arange(2, N + 1); pk = spf[2:]
m = n.copy()
while True:
    mask = (m % pk == 0) & (m > 1)
    if not mask.any(): break
    m = np.where(mask, m // pk, m)
Lam = np.zeros(N + 1); Lam[2:][m == 1] = np.log(pk[m == 1])
primes = n[pk == n].astype(float)
assert len(primes) == 78498
out = {"N": N, "primes": len(primes), "psi_N_minus_N": float(Lam.sum() - N)}

# zero heights (Odlyzko table), verified by the eta series
G = [14.134725141734693, 21.022039638771555, 25.010857580145688, 30.424876125859513, 32.935061587739189,
     37.586178158825671, 40.918719012147495, 43.327073280914999, 48.005150881167159, 49.773832477672302]
M = 400000
k = np.arange(1, M + 2, dtype=float); logk = np.log(k); sign = (-1.0) ** (k - 1)
def zeta_eta(s):
    S = np.cumsum(sign * np.exp(-s * logk)); eta = 0.5 * (S[-1] + S[-2]); return eta / (1 - 2 ** (1 - s))
zres = []
for g in G:
    z = zeta_eta(complex(0.5, g)); zres.append(abs(z))
print("zero heights verified via eta: max |zeta(1/2+i gamma_k)| =", f"{max(zres):.2e}")
assert max(zres) < 1e-6
out["zeros"] = [{"gamma": g, "abs_zeta": float(a)} for g, a in zip(G, zres)]
# off-line control: zeta at 0.5 + 15i is NOT small
out["control_abs_zeta_0.5+15i"] = float(abs(zeta_eta(complex(0.5, 15.0))))
print("control |zeta(1/2+15i)| =", f"{out['control_abs_zeta_0.5+15i']:.4f}")

def euler_partial(s, P):
    pp = primes[primes <= P]; return complex(np.prod(1.0 / (1.0 - pp ** (-s))))
CK = [10, 100, 1000, 10**4, 10**5, 3*10**5, 10**6]
out["s2"] = [{"P": P, "value": euler_partial(2.0, P).real, "error": euler_partial(2.0, P).real - math.pi**2/6} for P in CK]
def main_term(s, X): return (X ** (1 - s) / ((1 - s) * math.log(X))).real
out["line"] = {}
for gi in (0, 1):
    s = complex(0.5, G[gi]); rows = []
    for P in CK:
        v = euler_partial(s, P)
        rows.append({"P": P, "abs": abs(v), "arg_deg": math.degrees(np.angle(v)), "log_abs": math.log(abs(v)), "main": main_term(s, P)})
    out["line"][f"gamma{gi+1}"] = rows
    print(f"gamma_{gi+1}: |P_X| at checkpoints:", ", ".join(f"{r['abs']:.4f}" for r in rows))
coef = Lam[2:] - 1.0
out["compensated"] = {}
for sig in (0.9, 0.75, 0.6, 0.5, 0.4):
    part = np.cumsum(coef * n ** (-sig)) - 1.0   # n = 1 term is (0 - 1)
    out["compensated"][str(sig)] = [{"x": x, "sum": float(part[x - 2])} for x in (10**3, 10**4, 10**5, 10**6)]
    print(f"compensated s={sig}:", ", ".join(f"{part[x-2]:+.4f}" for x in (10**3, 10**4, 10**5, 10**6)))
json.dump(out, open("product_edge_receipts.json", "w"), indent=1)
print("wrote product_edge_receipts.json; ALL CHECKS PASSED")
