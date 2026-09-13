# Receipts — the geometry / rings / moats arc (2026-08-26 → 2026-09-13)

Every number quoted in the Angle & Length Atlas, Rooms VII–VIII, and the September threads comes from a
script in this folder. Each script asserts its own pre-stated expectations; a failing assert is a
finding, not a crash to hide. Corrections are kept on display in the scripts' docstrings and in the
Atlas §14 / thread ledgers rather than silently fixed.

**Standing fence (every file):** nothing here bears on the truth of the Riemann Hypothesis. Rigor rule 10,
the Davenport–Heilbronn gate: any candidate mechanism that survives with the D–H zero set substituted
for ζ's is a consequence of mirror symmetry alone and bears nothing. `positivity/forced_sign.py` computes
the D–H off-line zero 0.808517 + 85.699348i from scratch; `mobius/mobius_exponent_check.py` certifies
it by the argument principle (winding number 1 on radius 1e-3).

Run any file with `python3 <file>` (Python 3.9+, numpy; mpmath for the zero computations; scipy only for
`qubits/qutrit3_flips_checked.py`). Receipt JSONs are the outputs of the corresponding scripts.

| folder | file | what it establishes (grade) |
|---|---|---|
| (root) | verify_angles.py | asserts for every atlas angle: halving dial 90°/2^Ω(n), sopfr dial = 45° iff prime or 4, the three red rings (balance / bearing / Weil-leg) that ARE RH restated [T] |
| (root) | angle_dial.py, angle_dial.png | the one-page dial figure — Room VII's Python twin |
| (root) | reciprocal_check.py | w = 1/s sends the critical line to the circle \|w−1\| = 1 (the trace = norm circle); lattice points on it: 2, 1±i [T] |
| euler_product | euler_product_check.py | finite Euler products at ρ₁ swing 32× and do not converge; abscissa of Σ(Λ(n)−1)n^{−s} = Θ [T, classical]; the "never near 0" correction on display |
| euler_product | product_edge_receipts.py | the numbers on Room VIII's theorem card |
| moduli | mod_infinity.py, mod_units.py, mod_four.py | the modulus route: mod ∞−1 is no reduction, mod a unit is total collapse, mod 4 is the Gaussian class split [T] |
| moduli | eisenstein.py, eisenstein_fig.py, eisenstein.png | Z[ω]; L(s,χ₃) first zero 8.039737156 from scratch |
| moduli | gaussian_norm_count.py | Gaussian primes counted by norm vs li(x): 167, 1232, 9601, 78438 at 10³…10⁶ |
| positivity | trace_norm_rh.py | trace(w) = norm(w) ⟺ Re(1/w) = ½; pair form D = Σ(2β−1)²/(\|ρ\|²\|1−ρ̄\|²) ≥ 0, D = 0 ⟺ RH — a free square [T] |
| positivity | forced_sign.py | the Davenport–Heilbronn function from scratch; the off-line zero; the "forced sign" argument killed |
| positivity | proven_vs_needed.py, product_positivity.py | what the product proves (3-4-1 inequality → PNT, zero-free region) vs what is needed (Weil / Li positivity) |
| positivity | split_inert_halves.py | half-products over split / inert primes have branch points, not zeros: exponents (+½,+½) at ζ zeros, (+½,−½) at χ₄ zeros |
| sieve | prime_locator.py, gap_fill.py, two_receivers.py, gaps_from_below.py, gaussian_check_again.py (+ JSON) | Legendre exact (1,205 / 78,331), the 2e^{−γ} seam, Dickman/Buchstab two receivers, Jacobsthal 2…34, the 114-gap anatomy, the Gaussian lattice sieve 313,088 = 4·(prime ideals) + 4 |
| mobius | mobius_side.py (+ JSON) | Möbius chain to 1/ζ; M(10^k) = −1,1,2,−23,−48,212; Mertens/Pólya false on display; D–H has a Dirichlet inverse but no multiplicativity; NB constant = 2λ₁ + D |
| mobius | mobius_exponent_check.py (+ JSON) | 0.466 = 3-point LS slope, 0.4345 = pointwise sup (label collision, both right); D–H zero certified by winding number |
| moats | gaussian_moat.py, moat_rung2.py (+ png, JSONs) | Gaussian moats to step 4 replicated against Gethner–Wagon–Wick 1998 / Tsuchimura 2005 to 0.02%; largest prime-free disk radius 6.403 in the radius-2000 disk |
| zero_sums | zero_sum_dictionary.py (+ JSON) | rim angles → Li's λ_n (Keiper to 1e−9); trace/norm → the Nyman–Beurling constant 0.0461914; combination sums → exponent Θ; moats do not translate |
| qubits | two_qubit_houses.py | the prime 2's residue rings in Q(i)/Q(√5)/Q(√17) label one 2-qubit register: addition house-blind, multiplication sees the house (unit groups 2/3/1) [T] |
| qubits | eigenring_check.mts | audit of the Grok "Eigenring" app: multiplicative mode is anti-diffusion (sign flip); class labels only pick the graph |
| qubits | superstars.py (+ JSON) | smallest-factor usage shares with/without 2,3,5,7; prime-state Fourier peaks = Ramanujan sums \|μ(q)\|/φ(q), μ(4)=0 frequency reads the Gaussian census |
| qubits | qutrit3_flips_checked.py (+ JSON) | audit of the submitted qutrit script: mirrored-dof F1, ddof F5; nulls survive; side A is a size band |
| qubits | qutrit_triangle_mod7.py (+ JSON) | the true qutrit triangle is a cubic character (conductor 7); Chebyshev bias lives only in quadratic characters; zeros of the cubic L-function mod 7 |
| swarm_qc | swarm_qc.py, README.md, swarm_qc_receipts.json | the "multidimensional quantum computer" simulator: exact 4^N vs MPS bond budget, house gate sets, rebit encoding; P1, P2, P4, P5 pass, P3 refuted as stated |
| tools | publish_repo.py, publish_room8.py | the Pages publishers used for Rooms VII and VIII |
| handoffs | handoff_2026-09-11_angle-atlas-seam-moat.md | conversation-state transfer for the arc (supersedes earlier handoffs) |

Provenance: Hyperagent threads cmtux0xus04iw06addjeo409k (atlas, seam, moats) and cmtxc9ycb03zi07adhmd2g2rb
(adjudications, qubits, closure). Author of record: Travis Bergen; scripts written with the Riemann agent.
