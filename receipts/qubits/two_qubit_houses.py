#!/usr/bin/env python3
"""two_qubit_houses.py  (2026-09-12) — the i× rotation of Travis's "+ is house-blind, × sees the house".

The prime 2 in three quadratic houses:  ramified in Q(i)  [x^2+1 = (x+1)^2 mod 2],
inert in Q(sqrt5)  [x^2+x+1 irreducible mod 2],  split in Q(sqrt17)  [x^2+x = x(x+1) mod 2].
In each house the residue ring O_K/(2) has 4 elements -> it labels a 2-qubit register.
CLAIM (pre-stated): additive groups identical in all three houses (so the two-qubit Fourier
transform H⊗H is the same), multiplicative tables differ completely — unit-group orders 2, 3, 1.
Everything below is exhaustive enumeration over the 4 elements; asserts are the receipt.
"""
import itertools, numpy as np

# elements a + b*w as (a,b) with a,b in F_2; w^2 = P1*w + P0 (mod 2) per house
HOUSES = {
    "Q(i)      ramified  w=i,        w^2 = -1 = 1":        (0, 1),   # w^2 = 0*w + 1
    "Q(sqrt5)  inert     w=phi,      w^2 = w + 1":         (1, 1),   # w^2 = 1*w + 1
    "Q(sqrt17) split     w=(1+r17)/2, w^2 = w + 4 = w":    (1, 0),   # w^2 = 1*w + 0
}
E = [(a, b) for a in (0, 1) for b in (0, 1)]
idx = {e: i for i, e in enumerate(E)}

def add(x, y):
    return ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2)

def mul(x, y, P1, P0):
    a, b = x; c, d = y
    # (a + b w)(c + d w) = ac + (ad + bc) w + bd w^2,  w^2 = P1 w + P0
    const = (a * c + b * d * P0) % 2
    lin = (a * d + b * c + b * d * P1) % 2
    return (const, lin)

def perm_matrix(f):
    M = np.zeros((4, 4), dtype=int)
    for e in E:
        M[idx[f(e)], idx[e]] = 1
    return M

H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
HH = np.kron(H, H)

results = {}
add_tables = []
for name, (P1, P0) in HOUSES.items():
    add_tab = tuple(tuple(idx[add(x, y)] for y in E) for x in E)
    mul_tab = tuple(tuple(idx[mul(x, y, P1, P0)] for y in E) for x in E)
    one = (1, 0)
    units = [x for x in E if any(mul(x, y, P1, P0) == one for y in E)]
    nilp = [x for x in E if mul(x, x, P1, P0) == (0, 0)]
    idem = [x for x in E if mul(x, x, P1, P0) == x]
    assert all(add(x, x) == (0, 0) for x in E), "additive order is not 2 — additive group not (Z/2)^2"
    # additive gate set: translations t_e : x -> x + e  (these are the X-type gates on 2 qubits)
    add_gates = {e: perm_matrix(lambda x, e=e: add(x, e)) for e in E}
    # multiplicative gate set: m_u : x -> u*x for units u (permutations, hence unitary)
    mul_gates = {u: perm_matrix(lambda x, u=u: mul(u, x, P1, P0)) for u in units}
    # Fourier side: does H⊗H commute with the additive translations up to phase?  (diagonalises them)
    diag_ok = all(np.allclose(np.abs(HH @ G @ HH.T), np.eye(4)) for G in add_gates.values())
    add_tables.append(add_tab)
    results[name] = dict(units=units, n_units=len(units), nilpotents=nilp, idempotents=idem,
                         add_gates_are_pauli_X_group=all(np.array_equal(G @ G, np.eye(4, dtype=int)) for G in add_gates.values()),
                         HH_diagonalises_translations=diag_ok,
                         unit_gate_group_order=len(mul_gates), mul_table=mul_tab)

# receipts
assert len(set(add_tables)) == 1, "additive tables differ — claim FALSE"
mul_tabs = [r["mul_table"] for r in results.values()]
assert len(set(mul_tabs)) == 3, "some multiplication tables coincide — claim FALSE"
assert sorted(r["n_units"] for r in results.values()) == [1, 2, 3]

print("ADDITIVE: identical tables in all three houses; every nonzero element has order 2 => (Z/2)^2;")
print("          the four translations are {I, X⊗I, I⊗X, X⊗X}, and H⊗H diagonalises all of them (same Fourier side).")
for name, r in results.items():
    print(f"\n{name}")
    print(f"  units {r['units']}  (unit-gate group of order {r['unit_gate_group_order']})   nilpotents {r['nilpotents']}   idempotents {r['idempotents']}")
    print(f"  add gates are Pauli-X group: {r['add_gates_are_pauli_X_group']}   HH diagonalises translations: {r['HH_diagonalises_translations']}")
print("\nVERDICT: addition is house-blind (one table, one Fourier transform); multiplication sees the house")
print("         (three different tables; unit groups of order 2 = ramified, 3 = inert, 1 = split).")
