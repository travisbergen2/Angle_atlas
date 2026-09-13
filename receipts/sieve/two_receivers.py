"""two_receivers.py — receipts for Travis's 2026-09-10 21:05 observation:
"the composite-generating density from the small primes creates a structured
depletion, and there is a corresponding structure when you look at the
unresolved numbers from the other boundary."

Smooth side (generated from below, Dickman rho) vs rough side (unresolved
from above, Buchstab omega), both against the independence constant e^{-gamma}.
numpy only.
"""
import numpy as np, math, json, time
t0 = time.time()
EG = math.exp(-0.5772156649015329)          # e^{-gamma} = 0.561459...

# ---- Buchstab omega and Dickman rho by integrating their delay equations ----
h = 1e-4
U = np.arange(1.0, 8.0 + h, h)
om = np.empty_like(U); I = np.zeros_like(U)   # I(v) = int_1^v omega
for i, u in enumerate(U):
    if u <= 2.0 + 1e-12:
        om[i] = 1.0 / u
    else:
        j = int(round((u - 1.0 - 1.0) / h))    # index of u-1
        om[i] = (1.0 + I[j]) / u                # u*omega(u) = 1 + int_1^{u-1} omega
    if i > 0: I[i] = I[i-1] + 0.5 * h * (om[i] + om[i-1])
V = np.arange(0.0, 8.0 + h, h)
rho = np.ones_like(V); J = 0.0
for i, u in enumerate(V):
    if u <= 1.0 + 1e-12: rho[i] = 1.0; continue
    j = int(round((u - 1.0) / h))
    # rho(u) = 1 - int_1^u rho(t-1)/t dt  (trapezoid, sequential)
    J += 0.5 * h * (rho[j] / u + rho[j-1] / (u - h))
    rho[i] = 1.0 - J
def omega(u): return float(np.interp(u, U, om))
def dick(u):  return float(np.interp(u, V, rho))
R = {}
R["check_omega(3)=(1+ln2)/3"] = (round(omega(3), 6), round((1 + math.log(2)) / 3, 6))
R["check_rho(2)=1-ln2"] = (round(dick(2), 6), round(1 - math.log(2), 6))
R["e^-gamma"] = round(EG, 6)
R["omega_table"] = {str(u): {"omega": round(omega(u), 5), "omega*e^gamma": round(omega(u) / EG, 4), "rho": round(dick(u), 6)}
                    for u in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0)}
k = np.argmax(om[U > 2.0]); R["omega_first_max"] = (round(float(U[U > 2.0][k]), 3), round(float(om[U > 2.0][k]), 5), round(float(om[U > 2.0][k]) / EG, 4))
R["omega_min"] = (2.0, 0.5, round(0.5 / EG, 4))

# ---- sieve to 1e6: smooth (Psi) and rough (Phi) counts -----------------------
N = 10**6; n = np.arange(N + 1)
spf = np.zeros(N + 1, dtype=np.int64)
for p in range(2, int(N**0.5) + 1):
    if spf[p] == 0:
        idx = np.arange(p*p, N+1, p); m = spf[idx] == 0; spf[idx[m]] = p
isp = spf == 0; isp[:2] = False; spf[isp] = n[isp]; primes = np.nonzero(isp)[0]
gpf = np.zeros(N + 1, dtype=np.int64)
for p in primes: gpf[p::p] = p
pi = np.cumsum(isp)
rows = {}
for u in (1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0):
    y = N ** (1.0 / u)
    Phi = 1 + int(np.sum(spf[2:] > y))                 # unresolved at resolution y (incl. 1)
    Psi = 1 + int(np.sum(gpf[2:] <= y))                # generated from primes <= y (incl. 1)
    ind = N * float(np.prod(1 - 1.0 / primes[primes <= y].astype(float)))
    buch = N * omega(u) / math.log(y)
    nprimes_surv = int(pi[N] - pi[int(y)])
    rows[str(u)] = {"y": round(y, 1), "Phi_rough": Phi, "Buchstab_pred": round(buch), "independence_pred": round(ind),
                    "Phi/indep": round(Phi / ind, 4), "Phi/Buchstab": round(Phi / buch, 4),
                    "prime_fraction_of_unresolved": round(nprimes_surv / Phi, 4), "(1/u)/omega": round((1/u) / omega(u), 4),
                    "Psi_smooth": Psi, "Dickman_pred": round(N * dick(u)), "Psi/Dickman": round(Psi / (N * dick(u)), 4)}
R["x=1e6"] = rows
# CORRECTED PRE-STATEMENT (kept on display): the first draft asserted
# Phi + Psi - 1 == N at y = sqrt(x) ("rough + smooth cover everything once").
# FALSE: numbers with a factor <= y AND a factor > y are neither.  Count them.
mixed = N + 1 - rows["2.0"]["Phi_rough"] - rows["2.0"]["Psi_smooth"]
R["seam_u=2_mixed_count"] = mixed
R["seam_u=2_fractions"] = {"rough_only(=primes>y, +1)": round(rows["2.0"]["Phi_rough"] / N, 4),
                           "smooth_only": round(rows["2.0"]["Psi_smooth"] / N, 4),
                           "mixed(both_receivers)": round(mixed / N, 4)}
R["runtime_s"] = round(time.time() - t0, 1)
json.dump(R, open("two_receivers_receipts.json", "w"), indent=1)
for k_, v in R.items():
    if k_ == "x=1e6":
        for u, r in v.items(): print(u, r)
    else: print(k_, v)
assert abs(omega(3) - (1 + math.log(2)) / 3) < 1e-4
assert abs(dick(2) - (1 - math.log(2))) < 1e-4
assert mixed > 0 and rows["2.0"]["prime_fraction_of_unresolved"] == 1.0
print("ALL ASSERTS PASS")
