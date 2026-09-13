# swarm_qc — a "multidimensional quantum computer" simulator for the coupled-agent idea

**What it is.** A classical simulator of N agents, each holding a **dimension-4 site** (two qubits =
one qudit) whose basis is labeled by the residue ring of the prime 2 in a chosen quadratic house
(Gaussian Q(i) — ramified; golden Q(√5) — inert; Q(√17) — split). The global state lives in 4^N
dimensions. Two representations run side by side:

| representation | what it is | role |
|---|---|---|
| `Dense` | exact 4^N state vector | ground truth (N ≤ ~10–12 on a laptop) |
| `MPS` | matrix-product state, bond cap χ; agent *i* owns tensor A_i(α, s, β) | the **channel model**: the bond between neighbours is the only thing that carries entanglement, and log₂χ bounds the bits that can cross it |

**What it is not.** A quantum computer. Nothing here is hardware, nothing bears on the Riemann
Hypothesis, and no zeta data enters. The simulator's product is the *seam*: the depth at which a
bounded channel (χ) can no longer carry the global amplitude, measured against exact ground truth.

## Parameter map
| symbol | meaning | default |
|---|---|---|
| `n` | agents (sites); qubits = 2n | 8 |
| `chi` | bond dimension cap = channel budget; log₂χ bits per cut | {2, 8, 32, 256} |
| `depth` | brickwork layers | 10 |
| house | which residue ring labels the site: `gaussian`, `golden`, `split` | all three tested |
| `seed` | RNG seed for random/Clifford circuits | 7 |

Site index x = 2a + b ↔ ring element a + b·w (qubit 0 = a, qubit 1 = b).

## Gate set
- **Additive (house-blind):** translations x → x + e are the two-qubit Pauli-X group; the two-site
  coupling `CADD` |x⟩|y⟩ → |x⟩|x+y⟩ is ring addition (bitwise XOR) — identical in every house.
- **Multiplicative (house-sensitive):** unit multiplications y → u·y are permutations, hence unitary
  gates; the *generator* multiplication y → w·y is a unitary of order 2 in Q(i), order 3 in Q(√5),
  and in Q(√17) it is **idempotent of rank 2 — a projection, not a gate** (w is idempotent mod 2).
- **Universal extras:** H, S, T on either qubit of a site; Haar-random SU(2); `CZB` between sites; `SWAP`.
- **Rebit encoding:** ψ → [Re ψ; Im ψ] with a label qubit; U = A + iB → [[A, −B], [B, A]] (real
  orthogonal); multiplication by i → the 90° rotation J on the label qubit.

## Pre-stated expectations (in the module docstring, registered before the run)
P1 GHZ chain: 1 bit at every cut, χ = 2 exact. P2 random brickwork (N = 8): small χ fails within a
few layers, χ = 256 = 4⁴ matches exact at all depths. P3 Clifford circuit: same entanglement growth
as random (entanglement ≠ hardness; Clifford is stabilizer-simulable — Gottesman 1998,
Aaronson–Gottesman 2004, from memory). P4 generator gates as above. P5 rebit round trip < 1e-12.

## How to run
```
python3 swarm_qc.py          # self-tests (asserts) then experiments; receipts -> swarm_qc_receipts.json
```
Python 3 + numpy only. Runtime ≈ 1–3 minutes at the defaults (the χ = 256 arm dominates).
Change defaults by editing `experiments(n=8, depth=10, chis=(2, 8, 32, 256), seed=7)`.

## Design notes / honest limits
- MPS truncation is TEBD-style with alternating sweep direction so the orthogonality centre sits at
  the bond being cut; entropies are read after a full canonical sweep. Truncation is optimal only
  in canonical form, so reported χ requirements are upper bounds.
- Fidelity is always measured against the exact dense state; for N > ~12 that reference is gone and
  only truncation weight remains as a diagnostic.
- The residue-ring labeling puts the arithmetic where it can do work — in the multiplicative gate
  set — not in the coupling topology. Addition and the Fourier side (H⊗H) cannot see the house.
- Citations from memory (flag before any deposit): Vidal 2003 (efficient simulation of slightly
  entangled computation); Gottesman–Knill; Aaronson–Gottesman 2004; Rudolph–Grover 2002 and
  Aharonov 2003 (real-amplitude / rebit universality).
