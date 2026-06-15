"""
balance.py — Stoichiometric equation balancing.

Finds the smallest integer coefficients that balance a chemical equation by
solving M·x = 0 over the rationals, where M is the element-by-species matrix
(reactants positive, products negative) plus a row for net charge. Pure
Python (Fraction-based Gaussian elimination), no extra dependencies.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from functools import reduce

from rdkit import Chem

from .molecule import Molecule


def _lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b) if a and b else max(a, b)


def _species_vector(mol: Molecule) -> tuple[dict[str, int], int]:
    """Element counts (with H) and formal charge for a species."""
    return mol.atom_counts(), Chem.GetFormalCharge(mol.mol)


def _null_space_vector(matrix: list[list[Fraction]]) -> list[Fraction] | None:
    """Return one non-trivial solution x of matrix·x = 0, or None."""
    if not matrix or not matrix[0]:
        return None
    m = [row[:] for row in matrix]
    rows, cols = len(m), len(m[0])
    pivot_cols: list[int] = []
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        pivot_cols.append(c)
        r += 1
        if r == rows:
            break
    free_cols = [c for c in range(cols) if c not in pivot_cols]
    if not free_cols:
        return None
    free = free_cols[0]
    x = [Fraction(0)] * cols
    x[free] = Fraction(1)
    for ri, c in enumerate(pivot_cols):
        x[c] = -m[ri][free]
    return x


def balance(reactants: list[Molecule], products: list[Molecule]) -> dict | None:
    """Balance reactants -> products. Returns coefficients or None if impossible.

    Result: {"reactant_coeffs": [...], "product_coeffs": [...], "equation_fa": str}
    """
    species = reactants + products
    n_r = len(reactants)
    if not reactants or not products:
        return None

    vecs = [_species_vector(m) for m in species]
    elements = sorted({el for counts, _ in vecs for el in counts})

    rows: list[list[Fraction]] = []
    for el in elements:
        row = []
        for j, (counts, _) in enumerate(vecs):
            val = counts.get(el, 0)
            row.append(Fraction(val if j < n_r else -val))
        rows.append(row)
    # charge conservation row
    charge_row = [Fraction(c if j < n_r else -c) for j, (_, c) in enumerate(vecs)]
    if any(c != 0 for c in charge_row):
        rows.append(charge_row)

    sol = _null_space_vector(rows)
    if sol is None:
        return None

    # scale to integers
    denom = reduce(_lcm, (f.denominator for f in sol), 1)
    ints = [int(f * denom) for f in sol]
    g = reduce(gcd, (abs(i) for i in ints if i != 0), 0)
    if g == 0:
        return None
    ints = [i // g for i in ints]
    # normalise sign so the first reactant is positive
    if ints[0] < 0:
        ints = [-i for i in ints]
    # a valid balance has all coefficients strictly positive
    if any(i <= 0 for i in ints):
        return None

    r_coeffs = ints[:n_r]
    p_coeffs = ints[n_r:]

    def side(coeffs, mols):
        parts = []
        for c, m in zip(coeffs, mols):
            parts.append((f"{c} " if c != 1 else "") + m.formula)
        return " + ".join(parts)

    equation = f"{side(r_coeffs, reactants)} → {side(p_coeffs, products)}"
    return {
        "reactant_coeffs": r_coeffs,
        "product_coeffs": p_coeffs,
        "equation_fa": equation,
        "balanced": True,
    }
