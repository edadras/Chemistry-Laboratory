"""
redox.py — Oxidation states, electrochemistry & electrolysis.

Provides:
  * Oxidation-number assignment for every atom in a molecule.
  * Standard reduction potentials and galvanic-cell EMF (with spontaneity).
  * Electrolysis product prediction for common aqueous electrolytes
    (cathode / anode products).
"""
from __future__ import annotations

from rdkit import Chem

from .molecule import Molecule

# Standard reduction potentials E° (V) for common half-reactions.
STANDARD_POTENTIALS: dict[str, dict] = {
    "F2/F-":     {"e": 2.87, "fa": "F₂ + 2e⁻ → 2F⁻"},
    "Au3+/Au":   {"e": 1.50, "fa": "Au³⁺ + 3e⁻ → Au"},
    "Cl2/Cl-":   {"e": 1.36, "fa": "Cl₂ + 2e⁻ → 2Cl⁻"},
    "O2/H2O":    {"e": 1.23, "fa": "O₂ + 4H⁺ + 4e⁻ → 2H₂O"},
    "Ag+/Ag":    {"e": 0.80, "fa": "Ag⁺ + e⁻ → Ag"},
    "Fe3+/Fe2+": {"e": 0.77, "fa": "Fe³⁺ + e⁻ → Fe²⁺"},
    "Cu2+/Cu":   {"e": 0.34, "fa": "Cu²⁺ + 2e⁻ → Cu"},
    "H+/H2":     {"e": 0.00, "fa": "2H⁺ + 2e⁻ → H₂ (مرجع)"},
    "Pb2+/Pb":   {"e": -0.13, "fa": "Pb²⁺ + 2e⁻ → Pb"},
    "Ni2+/Ni":   {"e": -0.25, "fa": "Ni²⁺ + 2e⁻ → Ni"},
    "Fe2+/Fe":   {"e": -0.44, "fa": "Fe²⁺ + 2e⁻ → Fe"},
    "Zn2+/Zn":   {"e": -0.76, "fa": "Zn²⁺ + 2e⁻ → Zn"},
    "Al3+/Al":   {"e": -1.66, "fa": "Al³⁺ + 3e⁻ → Al"},
    "Mg2+/Mg":   {"e": -2.37, "fa": "Mg²⁺ + 2e⁻ → Mg"},
    "Na+/Na":    {"e": -2.71, "fa": "Na⁺ + e⁻ → Na"},
    "K+/K":      {"e": -2.93, "fa": "K⁺ + e⁻ → K"},
    "Li+/Li":    {"e": -3.04, "fa": "Li⁺ + e⁻ → Li"},
}

# Electronegativity (Pauling) for oxidation-number assignment.
_EN = {
    "F": 3.98, "O": 3.44, "Cl": 3.16, "N": 3.04, "Br": 2.96, "I": 2.66,
    "S": 2.58, "C": 2.55, "H": 2.20, "P": 2.19, "B": 2.04, "Si": 1.90,
    "Na": 0.93, "K": 0.82, "Mg": 1.31, "Ca": 1.00, "Fe": 1.83, "Cu": 1.90,
    "Zn": 1.65, "Al": 1.61,
}


def oxidation_states(mol: Molecule) -> dict:
    """Assign an oxidation number to every atom (electronegativity rule)."""
    m = Chem.AddHs(mol.mol)
    atoms = []
    for atom in m.GetAtoms():
        sym = atom.GetSymbol()
        ox = atom.GetFormalCharge()
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            osym = other.GetSymbol()
            if osym == sym:
                continue  # homonuclear bond contributes 0
            order = int(bond.GetBondTypeAsDouble()) or 1
            en_self = _EN.get(sym, 2.0)
            en_other = _EN.get(osym, 2.0)
            # more electronegative atom gains the electrons (−), other (+)
            ox += order * (1 if en_self < en_other else -1)
        atoms.append({"symbol": sym, "idx": atom.GetIdx(), "oxidation_state": ox})
    return {"atoms": atoms,
            "by_element": _summarise(atoms)}


def _summarise(atoms: list[dict]) -> dict:
    out: dict[str, list[int]] = {}
    for a in atoms:
        out.setdefault(a["symbol"], [])
        if a["oxidation_state"] not in out[a["symbol"]]:
            out[a["symbol"]].append(a["oxidation_state"])
    return {k: sorted(v) for k, v in out.items()}


def cell_potential(cathode: str, anode: str) -> dict:
    """E°cell = E°cathode − E°anode; positive ⇒ spontaneous galvanic cell."""
    if cathode not in STANDARD_POTENTIALS or anode not in STANDARD_POTENTIALS:
        return {"ok": False, "error_fa": "نیم‌واکنش ناشناخته است",
                "available": list(STANDARD_POTENTIALS.keys())}
    ec = STANDARD_POTENTIALS[cathode]["e"]
    ea = STANDARD_POTENTIALS[anode]["e"]
    emf = ec - ea
    spontaneous = emf > 0
    return {
        "ok": True,
        "cathode": cathode, "anode": anode,
        "e_cathode": ec, "e_anode": ea,
        "emf": round(emf, 2),
        "spontaneous": spontaneous,
        "summary_fa": (
            f"E°سل = {ec} − ({ea}) = {round(emf,2)} ولت — "
            f"{'خودبه‌خودی (سلول گالوانی)' if spontaneous else 'غیرخودبه‌خودی (نیازمند سلول الکترولیتی)'}"
        ),
    }


# Electrolysis of common aqueous electrolytes → (cathode product, anode product).
_ELECTROLYSIS = {
    "NaCl":  ("H₂ (آب احیا می‌شود، نه Na⁺)", "Cl₂ (یون کلرید اکسید می‌شود)"),
    "CuSO4": ("Cu (مس رسوب می‌کند)", "O₂ (آب اکسید می‌شود)"),
    "H2O":   ("H₂", "O₂"),
    "KI":    ("H₂", "I₂"),
    "AgNO3": ("Ag", "O₂"),
    "NaOH":  ("H₂", "O₂"),
    "CuCl2": ("Cu", "Cl₂"),
}


def electrolysis(electrolyte: str) -> dict:
    key = electrolyte.strip().replace(" ", "")
    # normalise a few common aliases
    aliases = {"کلریدسدیم": "NaCl", "نمک": "NaCl", "آب": "H2O", "water": "H2O",
               "سولفاتمس": "CuSO4"}
    key = aliases.get(key, key)
    if key not in _ELECTROLYSIS:
        return {"ok": False, "error_fa": f"الکترولیت {electrolyte} در پایگاه نیست",
                "available": list(_ELECTROLYSIS.keys())}
    cat, an = _ELECTROLYSIS[key]
    return {
        "ok": True, "electrolyte": key,
        "cathode_fa": cat, "anode_fa": an,
        "summary_fa": f"برق‌کافت محلول {key}: در کاتد ⟶ {cat}؛ در آند ⟶ {an}.",
    }
