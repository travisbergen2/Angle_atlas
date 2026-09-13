#!/usr/bin/env python3
"""
qutrit3_flips_checked.py — Travis's qutrit3_flips.py (2026-09-13), run as submitted EXCEPT for three
marked CHECK points that report BOTH the submitted computation and the corrected one:
  CHECK-1 (F1): the submitted test histograms theta_pi AND theta_pibar = pi/2 - theta_pi together; the
           6-bin histogram is then mirror-symmetric by construction (bin k == bin 5-k), the statistic is
           exactly 2x the honest 3-bin statistic and is compared to a chi^2 with 5 dof instead of 2.
           Corrected: theta_pi alone, 3 bins over (0, pi/4].
  CHECK-2 (F5): the docstring says dof 1 (f fitted) but scipy.stats.chisquare is called with ddof=0 (dof 2).
           Both p-values reported.
  CHECK-3 (F4 comment): the inline comment said P(+|split)=5/6, P(-|split)=1/6; the header and the code use
           2/3, 1/3 (correct under Hecke equidistribution on (0, pi/4]). Comment fixed.
Everything else is byte-for-byte the submitted logic.
"""
from __future__ import annotations
import json, math, time
import numpy as np
from scipy import stats

D = 3  # qutrit local dimension

GAMMA = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178,
         40.918719, 43.327073, 48.005151, 49.773832, 52.970321, 56.446248, 59.347044]

W = np.exp(2j * np.pi / 3)
F3 = (1 / np.sqrt(3)) * np.array([[1, 1, 1], [1, W, W ** 2], [1, W ** 2, W]])

def perm_mod3(f):
    P = np.zeros((9, 9))
    for x in range(3):
        for y in range(3):
            xp, yp = f(x, y); P[3 * xp + yp, 3 * x + y] = 1.0
    return P
CADD3 = perm_mod3(lambda x, y: (x, (x + y) % 3))

def haar_su3(rng):
    z = (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) / np.sqrt(2)
    q, r = np.linalg.qr(z); d = np.diag(r) / np.abs(np.diag(r))
    return q * d

class Dense:
    def __init__(self, n):
        self.n = n; self.psi = np.zeros(D ** n, complex); self.psi[0] = 1.0
    def apply_1(self, i, G):
        t = self.psi.reshape([D] * self.n)
        t = np.tensordot(G, t, axes=([1], [i])); t = np.moveaxis(t, 0, i); self.psi = t.reshape(-1)
    def apply_2(self, i, G):
        t = self.psi.reshape([D] * self.n); G4 = G.reshape(D, D, D, D)
        t = np.tensordot(G4, t, axes=([2, 3], [i, i + 1])); t = np.moveaxis(t, [0, 1], [i, i + 1])
        self.psi = t.reshape(-1)
    def cut_entropies(self):
        out = []
        for k in range(1, self.n):
            M = self.psi.reshape(D ** k, -1); s = np.linalg.svd(M, compute_uv=False)
            p = s ** 2; p = p[p > 1e-15]; out.append(float(-np.sum(p * np.log2(p))))
        return out

class MPS:
    def __init__(self, n, chi):
        self.n, self.chi = n, chi; self.A = []
        for _ in range(n):
            a = np.zeros((1, D, 1), complex); a[0, 0, 0] = 1.0; self.A.append(a)
        self.discarded = 0.0
    def apply_1(self, i, G): self.A[i] = np.einsum('ab,lbr->lar', G, self.A[i])
    def apply_2(self, i, G, absorb_right=True):
        Al, Ar = self.A[i], self.A[i + 1]; Dl, _, Dm = Al.shape; _, _, Dr = Ar.shape
        theta = np.einsum('lam,mbr->labr', Al, Ar).reshape(Dl, D * D, Dr)
        theta = np.einsum('xy,lyr->lxr', G, theta).reshape(Dl * D, D * Dr)
        U, S, Vh = np.linalg.svd(theta, full_matrices=False)
        k = max(1, min(self.chi, int(np.count_nonzero(S > 1e-13))))
        self.discarded += float(np.sum(S[k:] ** 2)); U, S, Vh = U[:, :k], S[:k], Vh[:k]
        if absorb_right:
            self.A[i] = U.reshape(Dl, D, k); self.A[i + 1] = (S[:, None] * Vh).reshape(k, D, Dr)
        else:
            self.A[i] = (U * S).reshape(Dl, D, k); self.A[i + 1] = Vh.reshape(k, D, Dr)
    def canonicalize(self):
        for i in range(self.n - 1):
            Dl, _, Dr = self.A[i].shape
            Q, R = np.linalg.qr(self.A[i].reshape(Dl * D, Dr))
            self.A[i] = Q.reshape(Dl, D, -1); self.A[i + 1] = np.einsum('ab,bcd->acd', R, self.A[i + 1])
        ents = [0.0] * (self.n - 1)
        for i in range(self.n - 1, 0, -1):
            Dl, _, Dr = self.A[i].shape
            U, S, Vh = np.linalg.svd(self.A[i].reshape(Dl, D * Dr), full_matrices=False)
            k = max(1, min(self.chi, int(np.count_nonzero(S > 1e-13))))
            self.discarded += float(np.sum(S[k:] ** 2)); U, S, Vh = U[:, :k], S[:k], Vh[:k]
            self.A[i] = Vh.reshape(k, D, Dr); self.A[i - 1] = np.einsum('lam,mk->lak', self.A[i - 1], U * S)
            p = S ** 2; p = p / p.sum(); p = p[p > 1e-15]; ents[i - 1] = float(-np.sum(p * np.log2(p)))
        return ents
    def to_dense(self):
        v = self.A[0]
        for i in range(1, self.n): v = np.einsum('l...m,mbr->l...br', v, self.A[i])
        return v.reshape(-1)
    def bond_dims(self): return [a.shape[2] for a in self.A[:-1]]

def fidelity(psi, phi):
    psi = psi / np.linalg.norm(psi); phi = phi / np.linalg.norm(phi)
    return float(abs(np.vdot(psi, phi)) ** 2)

def ghz3_circuit(n):
    layers = [[("1", 0, F3)]]
    layers += [[("2", i, CADD3)] for i in range(n - 1)]
    return layers
def random3_layers(n, depth, rng):
    layers = []
    for d in range(depth):
        L = [("1", i, haar_su3(rng)) for i in range(n)]
        bonds = range(n - 1) if d % 2 == 0 else range(n - 2, -1, -1)
        L += [("2", i, CADD3) for i in bonds]
        layers.append(L)
    return layers
def run_layers(sim, layers):
    for d, L in enumerate(layers):
        right = (d % 2 == 0)
        for kind, i, G in L:
            if kind == "1": sim.apply_1(i, G)
            else:
                if isinstance(sim, MPS): sim.apply_2(i, G, absorb_right=right)
                else: sim.apply_2(i, G)

def sieve_primes(N):
    is_p = np.ones(N + 1, bool); is_p[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if is_p[i]: is_p[i * i::i] = False
    return np.flatnonzero(is_p)

def two_squares(p):
    for b in range(1, int(math.isqrt(p)) + 1):
        a2 = p - b * b; a = math.isqrt(a2)
        if a * a == a2 and a >= b: return a, b
    raise ValueError(p)

def qc_cls(p):
    if p == 2: return "ram"
    return "split" if p % 4 == 1 else "inert"

def flip3(s): return (2 - s) % 3
def zsector_of(phase): return int(3 * phase / (2 * math.pi)) % 3

def build_number_theory(N, gamma1):
    primes = sieve_primes(N); rows = []
    for p in primes:
        p = int(p); c = qc_cls(p)
        phase = (gamma1 * math.log(p)) % (2 * math.pi)
        z = zsector_of(phase)
        row = dict(p=p, cls=c, zsector=z, phase=phase)
        if c == "split":
            a, b = two_squares(p); theta = math.atan2(b, a)
            row["theta_pi"] = theta; row["theta_pibar"] = math.pi / 2 - theta
            row["flips"] = [int(theta > math.pi / 6), int(theta > math.pi / 12)]
        elif c == "inert":
            row["flips"] = [int(p % 3 == 2)]
        else:
            row["flips"] = []
        row["n_flips"] = sum(row["flips"]); row["flip_parity"] = row["n_flips"] % 2
        zc = z
        for _ in range(row["flip_parity"]): zc = flip3(zc)
        row["z_composed"] = zc
        row["phase_signed"] = phase if row["flip_parity"] == 0 else (2 * math.pi - phase) % (2 * math.pi)
        rows.append(row)
    return rows

def run_tests():
    t0 = time.time(); R = {}
    assert np.allclose(F3.conj().T @ F3, np.eye(3)); assert np.allclose(CADD3.T @ CADD3, np.eye(9))
    for s in range(3):
        assert flip3(flip3(s)) == s
        ph = (s + 0.37) * (2 * math.pi) / 3
        assert zsector_of((2 * math.pi - ph) % (2 * math.pi)) == flip3(zsector_of(ph))
    rng = np.random.default_rng(1); n = 4; layers = random3_layers(n, 5, rng)
    Dd = Dense(n); M = MPS(n, chi=D ** (n // 2) * D); run_layers(Dd, layers); run_layers(M, layers)
    f = fidelity(Dd.psi, M.to_dense()); assert f > 1 - 1e-9, f
    ents_D = Dd.cut_entropies(); ents_M = M.canonicalize(); assert np.allclose(ents_D, ents_M, atol=1e-6)
    R["mps_vs_dense_untruncated_fidelity"] = f
    a, b = two_squares(5); assert (a, b) == (2, 1)
    R["two_squares_5"] = [a, b]; R["trace_5"] = 2 * a; R["norm_5"] = a * a + b * b
    R["tests_seconds"] = round(time.time() - t0, 2)
    return R

def experiments(n_qutrits=6, depth=8, chis=(3, 9, 27), seed=11, N_primes=200_000):
    R = {}
    lay = ghz3_circuit(n_qutrits); Dd = Dense(n_qutrits); run_layers(Dd, lay)
    M = MPS(n_qutrits, 3); run_layers(M, lay); ents = M.canonicalize()
    R["Q1_ghz3"] = dict(cut_entropies_bits=[round(e, 6) for e in ents], expected_bits=round(math.log2(3), 6),
                        chi=3, fidelity_vs_exact=fidelity(Dd.psi, M.to_dense()), max_bond=max(M.bond_dims()))
    rng = np.random.default_rng(seed); layers = random3_layers(n_qutrits, depth, rng)
    Dd = Dense(n_qutrits); dense_ents = []
    for L in layers:
        run_layers(Dd, [L]); dense_ents.append(Dd.cut_entropies())
    res = {"exact_cut_entropies_final": [round(e, 3) for e in dense_ents[-1]],
           "cap_bits_per_cut": [round(math.log2(3) * min(k, n_qutrits - k), 3) for k in range(1, n_qutrits)], "by_chi": {}}
    for chi in chis:
        Dc = Dense(n_qutrits); Mc = MPS(n_qutrits, chi); fids = []
        for d, L in enumerate(layers):
            run_layers(Dc, [L]); run_layers(Mc, [L])
            if d % 2 == 1 or d == len(layers) - 1: fids.append((d + 1, round(fidelity(Dc.psi, Mc.to_dense()), 6)))
        res["by_chi"][str(chi)] = dict(fidelity_by_depth=fids, first_depth_below_0p999=next((d for d, f in fids if f < 0.999), None),
                                      max_bond=max(Mc.bond_dims()))
    R["Q2_random3"] = res

    rows = build_number_theory(N_primes, GAMMA[0])
    split_rows = [r for r in rows if r["cls"] == "split"]; inert_rows = [r for r in rows if r["cls"] == "inert"]
    nonram_rows = split_rows + inert_rows
    for r in nonram_rows:
        ph = r["phase_signed"]
        if ph > 1e-9 and (2 * math.pi - ph) > 1e-9: assert zsector_of(ph) == r["z_composed"]

    # ---- F1  CHECK-1: submitted (mirrored, 6 bins over [0, pi/2)) vs corrected (theta_pi only, 3 bins over (0, pi/4])
    angles = []
    for r in split_rows: angles += [r["theta_pi"], r["theta_pibar"]]
    counts6, _ = np.histogram(angles, bins=np.linspace(0, math.pi / 2, 7)); chi2_6, p_6 = stats.chisquare(counts6)
    th = np.array([r["theta_pi"] for r in split_rows])
    counts3, _ = np.histogram(th, bins=np.linspace(0, math.pi / 4, 4)); chi2_3, p_3 = stats.chisquare(counts3)
    R["F1_hecke"] = dict(submitted_mirrored=dict(bin_counts=counts6.tolist(), chi2=round(float(chi2_6), 3), dof_assumed=5, p=round(float(p_6), 4)),
                         corrected_theta_only=dict(bin_counts=counts3.tolist(), chi2=round(float(chi2_3), 3), dof=2, p=round(float(p_3), 4),
                                                   reject_at_0p01=bool(p_3 < 0.01)),
                         note="mirrored statistic == 2 x theta-only statistic by construction: %.3f vs %.3f" % (chi2_6, 2 * chi2_3))

    classes = ["inert", "split"]; table = np.zeros((3, 2), int)
    for r in nonram_rows: table[r["zsector"], classes.index(r["cls"])] += 1
    c2, p2, dof2, _ = stats.chi2_contingency(table)
    R["F2_zetaphase_vs_class"] = dict(table=table.tolist(), chi2=round(float(c2), 3), dof=int(dof2), p=round(float(p2), 4), reject_at_0p01=bool(p2 < 0.01))
    table3 = np.zeros((3, 3), int)
    for r in nonram_rows: table3[r["zsector"], r["n_flips"]] += 1
    c3, p3, dof3, _ = stats.chi2_contingency(table3)
    R["F3_zetaphase_vs_flipcount"] = dict(table=table3.tolist(), chi2=round(float(c3), 3), dof=int(dof3), p=round(float(p3), 4), reject_at_0p01=bool(p3 < 0.01))

    # ---- F4 (CHECK-3: registered P(+|split)=2/3, P(-|split)=1/3 — the '5/6' comment in the submission was stale)
    def parity_counts(rr):
        plus = sum(1 for r in rr if r["flip_parity"] == 0); return plus, len(rr) - plus
    sp, sm = parity_counts(split_rows); ip, im = parity_counts(inert_rows)
    c_sp, p_sp = stats.chisquare([sp, sm], [len(split_rows) * 2 / 3, len(split_rows) / 3])
    c_ip, p_ip = stats.chisquare([ip, im], [len(inert_rows) / 2, len(inert_rows) / 2])
    R["F4_flip_alignment"] = dict(split=dict(plus=sp, minus=sm, expected_plus=round(len(split_rows) * 2 / 3, 1), p=round(float(p_sp), 4)),
                                  inert=dict(plus=ip, minus=im, expected_plus=round(len(inert_rows) / 2, 1), p=round(float(p_ip), 4)))

    # ---- F5  CHECK-2: ddof=0 (as submitted) vs ddof=1 (as the docstring states)
    raw_hist = [sum(1 for r in nonram_rows if r["zsector"] == s) for s in range(3)]
    comp_hist = [sum(1 for r in nonram_rows if r["z_composed"] == s) for s in range(3)]
    f_odd = sum(1 for r in nonram_rows if r["flip_parity"] == 1) / len(nonram_rows)
    exp_comp = [(1 - f_odd) * raw_hist[0] + f_odd * raw_hist[2], float(raw_hist[1]), (1 - f_odd) * raw_hist[2] + f_odd * raw_hist[0]]
    c5, p5_ddof0 = stats.chisquare(comp_hist, exp_comp); _, p5_ddof1 = stats.chisquare(comp_hist, exp_comp, ddof=1)
    R["F5_composed_vs_mirror_mix"] = dict(raw_sector_hist=raw_hist, composed_sector_hist=comp_hist, flip_rate_odd=round(f_odd, 4),
                                          expected=[round(e, 1) for e in exp_comp], chi2=round(float(c5), 3),
                                          p_submitted_ddof0=round(float(p5_ddof0), 4), p_stated_dof1=round(float(p5_ddof1), 4))
    # ---- extra diagnostic (not in the submission): is the 'zeta' sector uniform? Weyl sum |sum p^{i gamma}|/pi(x) vs 1/sqrt(1+gamma^2)
    ps = np.array([r["p"] for r in rows], float)
    weyl = abs(np.exp(1j * GAMMA[0] * np.log(ps)).sum()) / len(ps)
    R["diag_sideA"] = dict(raw_sector_hist_all_primes=[sum(1 for r in rows if r["zsector"] == s) for s in range(3)],
                           uniform_p=round(float(stats.chisquare([sum(1 for r in rows if r["zsector"] == s) for s in range(3)])[1]), 6),
                           weyl_sum_over_pi=round(float(weyl), 4), predicted_1_over_sqrt_1_plus_gamma2=round(1 / math.sqrt(1 + GAMMA[0] ** 2), 4),
                           sector_band_ratio_in_p=round(math.exp(2 * math.pi / (3 * GAMMA[0])), 4))
    R["n_split"], R["n_inert"] = len(split_rows), len(inert_rows)
    return R

if __name__ == "__main__":
    out = {"tests": run_tests()}; print("TESTS PASSED:", json.dumps(out["tests"]))
    t0 = time.time(); ex = out["experiments"] = experiments(); out["seconds"] = round(time.time() - t0, 1)
    print("Q1", ex["Q1_ghz3"])
    print("Q2 cap", ex["Q2_random3"]["cap_bits_per_cut"], "exact final", ex["Q2_random3"]["exact_cut_entropies_final"])
    for chi, v in ex["Q2_random3"]["by_chi"].items(): print("   chi", chi, v)
    for k in ("F1_hecke", "F2_zetaphase_vs_class", "F3_zetaphase_vs_flipcount", "F4_flip_alignment", "F5_composed_vs_mirror_mix", "diag_sideA"):
        print(k, json.dumps(ex[k]))
    print("n_split", ex["n_split"], "n_inert", ex["n_inert"], "seconds", out["seconds"])
    json.dump(out, open("qutrit3_flips_checked_receipts.json", "w"), indent=1)
