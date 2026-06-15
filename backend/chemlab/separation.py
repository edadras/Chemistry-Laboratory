"""
separation.py — Distillation & separation simulation.

Simulates separating a liquid mixture by boiling point: ranks components,
decides which can be separated by simple vs. fractional distillation (based on
the boiling-point gap), and reports the order in which fractions come over.
Boiling points come from a curated table; unknowns are estimated from molar
mass and polarity as a rough fallback.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from .molecule import Molecule

# Boiling points (°C, 1 atm) for common liquids, keyed by canonical SMILES.
BOILING_POINT_C: dict[str, float] = {
    "O": 100.0,            # water
    "CO": 64.7,            # methanol
    "CCO": 78.37,          # ethanol
    "CC(C)O": 82.6,        # isopropanol
    "CC(=O)C": 56.0,       # acetone
    "CC(=O)O": 118.0,      # acetic acid
    "CCOC(C)=O": 77.1,     # ethyl acetate
    "c1ccccc1": 80.1,      # benzene
    "Cc1ccccc1": 110.6,    # toluene
    "ClCCl": 39.6,         # DCM-ish (dichloromethane is ClCCl? actually C(Cl)Cl)
    "CCCCCC": 69.0,        # hexane
    "CCCCC": 36.1,         # pentane
    "ClC(Cl)Cl": 61.2,     # chloroform
    "CCOCC": 34.6,         # diethyl ether
    "OCC(O)CO": 290.0,     # glycerol
    "CC#N": 82.0,          # acetonitrile
    "CS(C)=O": 189.0,      # DMSO
}


def boiling_point(mol: Molecule) -> tuple[float, bool]:
    """Return (bp_celsius, is_estimated)."""
    try:
        canon = Chem.CanonSmiles(mol.smiles)
    except Exception:
        canon = mol.smiles
    if canon in BOILING_POINT_C:
        return BOILING_POINT_C[canon], False
    # rough estimate: heavier & more polar → higher bp
    mw = Descriptors.MolWt(mol.mol)
    logp = Crippen.MolLogP(mol.mol)
    est = -50 + 1.6 * mw - 8 * logp
    return round(est, 1), True


def distill(components: list[Molecule]) -> dict:
    """Plan a distillation separation of a mixture by boiling point."""
    if len(components) < 1:
        return {"ok": False, "error_fa": "حداقل یک جزء لازم است"}
    rows = []
    for m in components:
        bp, est = boiling_point(m)
        rows.append({"formula": m.formula, "smiles": m.smiles,
                     "bp_c": bp, "estimated": est})
    rows.sort(key=lambda r: r["bp_c"])

    # decide separability from boiling-point gaps
    notes = []
    method = "تقطیر ساده" if len(rows) > 1 else "بدون نیاز به جداسازی (تک‌جزء)"
    min_gap = None
    for a, b in zip(rows, rows[1:]):
        gap = b["bp_c"] - a["bp_c"]
        min_gap = gap if min_gap is None else min(min_gap, gap)
    if min_gap is not None:
        if min_gap < 25:
            method = "تقطیر جزءبه‌جزء (اختلاف نقطه‌ی جوش کم)"
            notes.append(f"کمترین اختلاف نقطه‌ی جوش ≈ {round(min_gap,1)}°C — نیاز به ستون تقطیر جزءبه‌جزء")
        else:
            notes.append(f"کمترین اختلاف نقطه‌ی جوش ≈ {round(min_gap,1)}°C — تقطیر ساده کافی است")

    order = "، سپس ".join(
        f"{r['formula']} (~{r['bp_c']}°C)" for r in rows)
    return {
        "ok": True,
        "method_fa": method,
        "fractions": rows,  # in order they distil over (low bp first)
        "min_gap_c": round(min_gap, 1) if min_gap is not None else None,
        "notes_fa": notes,
        "summary_fa": f"روش پیشنهادی: {method}. ترتیب خروج تقطیر (از سبک به سنگین): {order}.",
    }
