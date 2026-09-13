#!/usr/bin/env python3
"""
swarm_qc.py — "multidimensional quantum computer" SIMULATOR for the coupled-agent idea (2026-09-12).

N agents, each holding a dimension-4 site (two qubits = one qudit) whose basis is labeled by the
residue ring of the prime 2 in a chosen quadratic house:
    gaussian  Q(i)      2 ramified   w = i,   w^2 = 1 (mod 2)          ring F2[x]/(x+1)^2
    golden    Q(sqrt5)  2 inert      w = phi, w^2 = w + 1              field F4
    split     Q(sqrt17) 2 split      w = (1+sqrt17)/2, w^2 = w        ring F2 x F2
Element a + b*w  <->  site basis index 2a + b  (qubit 0 = a, qubit 1 = b).

Two representations of the global 4^N-dimensional state:
    Dense  — exact state vector (ground truth; feasible to N ~ 10-12 here)
    MPS    — matrix-product state with bond cap chi: each agent owns a tensor A_i(alpha, s, beta);
             the bond is the channel between neighbouring agents; log2(chi) bounds the entanglement
             (in bits) that can cross it.  Two-site gates: contract -> apply -> SVD -> truncate.
             Layers are applied as alternating left/right sweeps so the orthogonality centre sits
             at the bond being cut (TEBD-style); entropies are read after a full canonical sweep.

PRE-STATED EXPECTATIONS (registered before running; beatable):
  P1  GHZ chain across all sites: every cut carries exactly 1 bit; chi = 2 reproduces the exact
      state to < 1e-12.                                                   [theorem: Schmidt rank 2]
  P2  Random brickwork, N = 8: cut entropies climb toward the cap min(k, N-k)*2 bits; chi = 2 and
      chi = 8 drop below fidelity 0.999 within a few layers; chi = 256 = 4^4 (exact middle-cut
      maximum) matches the dense state to ~1e-10 at every depth.           [prediction]
  P3  Clifford-only brickwork shows the same entropy growth as the random circuit within ~25%
      at the final depth — entanglement alone does not predict hardness (Clifford circuits are
      efficiently simulable by stabilizer methods: Gottesman 1998 / Aaronson-Gottesman 2004,
      from memory — flag; not implemented here).                          [prediction, weak]
  P4  Generator multiplication M_w: y -> w*y is a permutation (unitary) of order 2 in Q(i) and
      order 3 in Q(sqrt5); in Q(sqrt17) it is idempotent of rank 2 — a projection, not a gate.
                                                                          [theorem by enumeration]
  P5  Rebit encoding: psi -> [Re psi; Im psi] with a label qubit; U = A + iB -> [[A,-B],[B,A]]
      (real orthogonal); i -> the 90-degree rotation J on the label qubit; round trip < 1e-12.
                                                                          [theorem]
Fence: this is classical simulation of quantum circuits. Nothing here is a quantum computer,
nothing bears on RH, and no zeta data enters. Run:  python3 swarm_qc.py   (tests, then experiments;
receipts to swarm_qc_receipts.json).
"""
from __future__ import annotations
import json, sys, time
import numpy as np

# ----------------------------------------------------------------------------- houses
HOUSES = {"gaussian": (0, 1), "golden": (1, 1), "split": (1, 0)}   # w^2 = P1*w + P0 (mod 2)
E = [(0, 0), (0, 1), (1, 0), (1, 1)]                                # (a, b) <-> a + b w
def idx(e): return 2 * e[0] + e[1]
def ring_add(x, y): return ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2)
def ring_mul(x, y, house):
    P1, P0 = HOUSES[house]; a, b = x; c, d = y
    return ((a * c + b * d * P0) % 2, (a * d + b * c + b * d * P1) % 2)
def units(house):
    return [x for x in E if any(ring_mul(x, y, house) == (1, 0) for y in E)]
def gen_mul_matrix(house):
    """M_w : |y> -> |w*y| as a 4x4 real matrix (w = (0,1))."""
    M = np.zeros((4, 4))
    for y in E:
        M[idx(ring_mul((0, 1), y, house)), idx(y)] = 1.0
    return M
def unit_mul_matrix(u, house):
    M = np.zeros((4, 4))
    for y in E:
        M[idx(ring_mul(u, y, house)), idx(y)] = 1.0
    return M

# ----------------------------------------------------------------------------- gates
I2 = np.eye(2); I4 = np.eye(4)
H1 = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S1 = np.diag([1, 1j]); T1 = np.diag([1, np.exp(1j * np.pi / 4)])
X1 = np.array([[0, 1], [1, 0]])
def on_q0(g): return np.kron(g, I2)          # qubit 0 = high bit of the site index
def on_q1(g): return np.kron(I2, g)
def haar_su2(rng):
    z = (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))) / np.sqrt(2)
    q, r = np.linalg.qr(z); d = np.diag(r) / np.abs(np.diag(r))
    return q * d
def perm16(f):
    """16x16 permutation from a map (x, y) -> (x', y') on site-index pairs."""
    P = np.zeros((16, 16))
    for x in range(4):
        for y in range(4):
            xp, yp = f(x, y); P[4 * xp + yp, 4 * x + y] = 1.0
    return P
CADD = perm16(lambda x, y: (x, x ^ y))            # |x>|y> -> |x>|x+y>  (ring addition = XOR; house-blind)
SWAP = perm16(lambda x, y: (y, x))
def cz_between(): # CZ between qubit 1 of the left site and qubit 0 of the right site
    d = np.ones(16)
    for x in range(4):
        for y in range(4):
            if (x & 1) and (y >> 1): d[4 * x + y] = -1
    return np.diag(d)
CZB = cz_between()

# ----------------------------------------------------------------------------- dense reference
class Dense:
    def __init__(self, n):
        self.n = n; self.psi = np.zeros(4 ** n, complex); self.psi[0] = 1.0
    def apply_1(self, i, G):
        t = self.psi.reshape([4] * self.n)
        t = np.tensordot(G, t, axes=([1], [i])); t = np.moveaxis(t, 0, i)
        self.psi = t.reshape(-1)
    def apply_2(self, i, G):
        t = self.psi.reshape([4] * self.n); G4 = G.reshape(4, 4, 4, 4)
        t = np.tensordot(G4, t, axes=([2, 3], [i, i + 1])); t = np.moveaxis(t, [0, 1], [i, i + 1])
        self.psi = t.reshape(-1)
    def cut_entropies(self):
        out = []
        for k in range(1, self.n):
            M = self.psi.reshape(4 ** k, -1); s = np.linalg.svd(M, compute_uv=False)
            p = s ** 2; p = p[p > 1e-15]; out.append(float(-np.sum(p * np.log2(p))))
        return out

# ----------------------------------------------------------------------------- MPS
class MPS:
    def __init__(self, n, chi):
        self.n, self.chi = n, chi
        self.A = []
        for _ in range(n):
            a = np.zeros((1, 4, 1), complex); a[0, 0, 0] = 1.0; self.A.append(a)
        self.discarded = 0.0
    def apply_1(self, i, G):
        self.A[i] = np.einsum('ab,lbr->lar', G, self.A[i])
    def apply_2(self, i, G, absorb_right=True):
        Al, Ar = self.A[i], self.A[i + 1]
        Dl, _, Dm = Al.shape; _, _, Dr = Ar.shape
        theta = np.einsum('lam,mbr->labr', Al, Ar).reshape(Dl, 16, Dr)
        theta = np.einsum('xy,lyr->lxr', G, theta).reshape(Dl * 4, 4 * Dr)
        U, S, Vh = np.linalg.svd(theta, full_matrices=False)
        k = max(1, min(self.chi, int(np.count_nonzero(S > 1e-13))))
        self.discarded += float(np.sum(S[k:] ** 2)); U, S, Vh = U[:, :k], S[:k], Vh[:k]
        if absorb_right:
            self.A[i] = U.reshape(Dl, 4, k); self.A[i + 1] = (S[:, None] * Vh).reshape(k, 4, Dr)
        else:
            self.A[i] = (U * S).reshape(Dl, 4, k); self.A[i + 1] = Vh.reshape(k, 4, Dr)
    def canonicalize(self):
        """Left QR sweep then right SVD sweep; returns Schmidt entropies (bits) per bond."""
        for i in range(self.n - 1):
            Dl, _, Dr = self.A[i].shape
            Q, R = np.linalg.qr(self.A[i].reshape(Dl * 4, Dr))
            self.A[i] = Q.reshape(Dl, 4, -1); self.A[i + 1] = np.einsum('ab,bcd->acd', R, self.A[i + 1])
        ents = [0.0] * (self.n - 1)
        for i in range(self.n - 1, 0, -1):
            Dl, _, Dr = self.A[i].shape
            U, S, Vh = np.linalg.svd(self.A[i].reshape(Dl, 4 * Dr), full_matrices=False)
            k = max(1, min(self.chi, int(np.count_nonzero(S > 1e-13))))
            self.discarded += float(np.sum(S[k:] ** 2)); U, S, Vh = U[:, :k], S[:k], Vh[:k]
            self.A[i] = Vh.reshape(k, 4, Dr); self.A[i - 1] = np.einsum('lam,mk->lak', self.A[i - 1], U * S)
            p = S ** 2; p = p / p.sum(); p = p[p > 1e-15]; ents[i - 1] = float(-np.sum(p * np.log2(p)))
        return ents
    def to_dense(self):
        v = self.A[0]
        for i in range(1, self.n):
            v = np.einsum('l...m,mbr->l...br', v, self.A[i])
        return v.reshape(-1)
    def bond_dims(self): return [a.shape[2] for a in self.A[:-1]]

def fidelity(psi, phi):
    psi = psi / np.linalg.norm(psi); phi = phi / np.linalg.norm(phi)
    return float(abs(np.vdot(psi, phi)) ** 2)

# ----------------------------------------------------------------------------- circuits (layer lists)
def ghz_circuit(n):
    layers = [[("1", 0, on_q0(H1))]]
    layers += [[("2", i, CADD)] for i in range(n - 1)]      # |x>|0> -> |x>|x> down the chain
    return layers
def random_layers(n, depth, rng):
    layers = []
    for d in range(depth):
        L = [("1", i, on_q0(haar_su2(rng)) @ on_q1(haar_su2(rng))) for i in range(n)]
        bonds = range(n - 1) if d % 2 == 0 else range(n - 2, -1, -1)
        L += [("2", i, CADD if rng.random() < 0.5 else CZB) for i in bonds]
        layers.append(L)
    return layers
def clifford_layers(n, depth, rng):
    cliff1 = [on_q0(H1), on_q1(H1), on_q0(S1), on_q1(S1), on_q0(H1) @ on_q1(S1), on_q1(H1) @ on_q0(S1)]
    layers = []
    for d in range(depth):
        L = [("1", i, cliff1[rng.integers(len(cliff1))]) for i in range(n)]
        bonds = range(n - 1) if d % 2 == 0 else range(n - 2, -1, -1)
        L += [("2", i, CADD if rng.random() < 0.5 else CZB) for i in bonds]
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

# ----------------------------------------------------------------------------- rebit encoding
def rebit_encode(psi): return np.concatenate([psi.real, psi.imag])
def rebit_decode(r): D = r.size // 2; return r[:D] + 1j * r[D:]
def rebit_gate(U):
    A, B = U.real, U.imag
    return np.block([[A, -B], [B, A]])
def rebit_i(D):  # multiplication by i on C^D = rotation J on the label qubit
    return np.kron(np.array([[0.0, -1.0], [1.0, 0.0]]), np.eye(D))

# ----------------------------------------------------------------------------- tests
def run_tests():
    rng = np.random.default_rng(1)
    t0 = time.time(); R = {}
    # gates unitary / permutations
    for G in (CADD, SWAP, CZB, on_q0(H1), on_q1(S1)):
        assert np.allclose(G.conj().T @ G, np.eye(G.shape[0]))
    # dense vs MPS with no truncation, random circuit n=5
    n = 5; layers = random_layers(n, 6, rng)
    D = Dense(n); M = MPS(n, chi=4 ** (n // 2) * 4); run_layers(D, layers); run_layers(M, layers)
    assert abs(np.linalg.norm(D.psi) - 1) < 1e-10
    f = fidelity(D.psi, M.to_dense()); assert f > 1 - 1e-9, f
    ents_D = D.cut_entropies(); ents_M = M.canonicalize()
    assert np.allclose(ents_D, ents_M, atol=1e-6), (ents_D, ents_M)
    R["mps_vs_dense_untruncated_fidelity"] = f
    # P4 house generators
    R["P4"] = {}
    for h in HOUSES:
        Mw = gen_mul_matrix(h); sv = np.linalg.svd(Mw, compute_uv=False)
        unitary = np.allclose(Mw.T @ Mw, I4)
        order = None
        if unitary:
            P = np.eye(4);
            for k in range(1, 7):
                P = P @ Mw
                if np.allclose(P, np.eye(4)): order = k; break
        R["P4"][h] = dict(unitary=bool(unitary), order=order, idempotent=bool(np.allclose(Mw @ Mw, Mw)),
                          rank=int(np.linalg.matrix_rank(Mw)), singular_values=[round(float(s), 6) for s in sv],
                          n_units=len(units(h)))
    assert R["P4"]["gaussian"]["unitary"] and R["P4"]["gaussian"]["order"] == 2
    assert R["P4"]["golden"]["unitary"] and R["P4"]["golden"]["order"] == 3
    assert (not R["P4"]["split"]["unitary"]) and R["P4"]["split"]["idempotent"] and R["P4"]["split"]["rank"] == 2
    # P5 rebit
    psi = rng.normal(size=16) + 1j * rng.normal(size=16); psi /= np.linalg.norm(psi)
    U = np.kron(np.kron(haar_su2(rng), haar_su2(rng)), np.kron(haar_su2(rng), haar_su2(rng)))  # 16x16 unitary
    UR = rebit_gate(U); assert np.allclose(UR.T @ UR, np.eye(32))
    err = np.linalg.norm(rebit_decode(UR @ rebit_encode(psi)) - U @ psi)
    err_i = np.linalg.norm(rebit_decode(rebit_i(16) @ rebit_encode(psi)) - 1j * psi)
    assert err < 1e-12 and err_i < 1e-12 and abs(np.linalg.norm(rebit_encode(psi)) - 1) < 1e-12
    R["P5"] = dict(roundtrip_error=float(err), i_as_rotation_error=float(err_i), UR_orthogonal=True)
    # unit gates are permutations in every house
    for h in HOUSES:
        for u in units(h):
            Mu = unit_mul_matrix(u, h); assert np.allclose(Mu.T @ Mu, I4)
    R["tests_seconds"] = round(time.time() - t0, 2)
    return R

# ----------------------------------------------------------------------------- experiments
def experiments(n=8, depth=10, chis=(2, 8, 32, 256), seed=7):
    rng = np.random.default_rng(seed); R = {"n_sites": n, "n_qubits": 2 * n, "depth": depth}
    # P1 GHZ
    lay = ghz_circuit(n); D = Dense(n); run_layers(D, lay)
    M = MPS(n, 2); run_layers(M, lay); ents = M.canonicalize()
    R["P1_ghz"] = dict(cut_entropies_bits=[round(e, 6) for e in ents], chi=2,
                       fidelity_vs_exact=fidelity(D.psi, M.to_dense()), max_bond=max(M.bond_dims()))
    # P2 random, P3 clifford
    for name, maker in (("P2_random", random_layers), ("P3_clifford", clifford_layers)):
        layers = maker(n, depth, np.random.default_rng(seed + (0 if name.startswith("P2") else 1)))
        D = Dense(n); dense_ents = []
        for L in layers:
            run_layers(D, [L]); dense_ents.append(D.cut_entropies())
        # rerun dense layer by layer for fidelity checkpoints per chi
        res = {"exact_cut_entropies_by_depth": [[round(e, 3) for e in row] for row in dense_ents],
               "cap_bits_per_cut": [2 * min(k, n - k) for k in range(1, n)], "by_chi": {}}
        for chi in chis:
            Dc = Dense(n); Mc = MPS(n, chi); fids = []; t0 = time.time()
            for d, L in enumerate(layers):
                run_layers(Dc, [L]); run_layers(Mc, [L])
                if d % 2 == 1 or d == len(layers) - 1:
                    fids.append((d + 1, round(fidelity(Dc.psi, Mc.to_dense()), 6)))
            first_fail = next((d for d, f in fids if f < 0.999), None)
            res["by_chi"][str(chi)] = dict(fidelity_by_depth=fids, first_depth_below_0p999=first_fail,
                                          discarded_weight=round(Mc.discarded, 6), max_bond=max(Mc.bond_dims()),
                                          seconds=round(time.time() - t0, 1))
        R[name] = res
    # P3 comparison: middle-cut entropy at final depth
    mid = n // 2 - 1
    e_rand = R["P2_random"]["exact_cut_entropies_by_depth"][-1][mid]
    e_cliff = R["P3_clifford"]["exact_cut_entropies_by_depth"][-1][mid]
    R["P3_vs_P2_middle_cut_final_bits"] = dict(random=e_rand, clifford=e_cliff, cap=2 * (n // 2),
                                                ratio=round(e_cliff / e_rand, 3) if e_rand else None)
    return R

if __name__ == "__main__":
    out = {"tests": run_tests()}
    print("TESTS PASSED:", json.dumps(out["tests"], indent=None)[:600])
    t0 = time.time(); out["experiments"] = experiments(); out["experiments_seconds"] = round(time.time() - t0, 1)
    ex = out["experiments"]
    print(f"\nP1 GHZ (N={ex['n_sites']} sites = {ex['n_qubits']} qubits): cut entropies {ex['P1_ghz']['cut_entropies_bits']} bits;"
          f" chi=2 fidelity {ex['P1_ghz']['fidelity_vs_exact']:.12f}")
    for name in ("P2_random", "P3_clifford"):
        r = ex[name]; print(f"\n{name}: cap bits per cut {r['cap_bits_per_cut']}")
        print("  exact cut entropies at final depth:", r["exact_cut_entropies_by_depth"][-1])
        for chi, v in r["by_chi"].items():
            print(f"  chi={chi:>3}: fidelity by depth {v['fidelity_by_depth']}  first<0.999 at depth {v['first_depth_below_0p999']}  max bond {v['max_bond']}  ({v['seconds']}s)")
    print("\nP3 vs P2 middle cut, final depth:", ex["P3_vs_P2_middle_cut_final_bits"])
    print("P4 generators:", json.dumps(out["tests"]["P4"]))
    print(f"\nexperiments {out['experiments_seconds']}s")
    json.dump(out, open("swarm_qc_receipts.json", "w"), indent=1)
