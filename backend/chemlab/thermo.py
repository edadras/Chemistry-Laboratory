"""
thermo.py — Reaction thermodynamics & kinetics (estimates).

Thermodynamics
--------------
ΔH is estimated with the *average bond-enthalpy* method:
    ΔH ≈ Σ E(bonds in reactants) − Σ E(bonds in products)
weighted by stoichiometric coefficients. Bonds unchanged on both sides cancel,
so spectator/aromatic inaccuracies wash out. Entropy is a crude sign estimate
from the change in molecule count, giving an approximate ΔG and a spontaneity
verdict. Values are estimates for teaching, not calorimetric data.

Kinetics
--------
A transparent rule-based estimate of reaction rate and yield from the reaction
type and the conditions (temperature via a qualitative Arrhenius factor,
catalyst, reversibility).
"""
from __future__ import annotations

import math

from rdkit import Chem

from .molecule import Molecule
from .conditions import Conditions

R = 8.314  # J/mol·K

# Average bond enthalpies (kJ/mol). Key: (symbol1, symbol2 sorted, order),
# order is 1/2/3 for single/double/triple and "ar" for aromatic.
BOND_ENTHALPY: dict[tuple, float] = {
    ("H", "H", 1): 436,
    ("C", "H", 1): 413, ("C", "C", 1): 348, ("C", "C", 2): 614, ("C", "C", 3): 839,
    ("C", "O", 1): 358, ("C", "O", 2): 799,
    ("H", "O", 1): 463, ("O", "O", 1): 146, ("O", "O", 2): 498,
    ("C", "N", 1): 308, ("C", "N", 2): 615, ("C", "N", 3): 891,
    ("H", "N", 1): 391, ("N", "N", 1): 163, ("N", "N", 2): 418, ("N", "N", 3): 941,
    ("N", "O", 1): 201, ("N", "O", 2): 607,
    ("C", "Cl", 1): 328, ("Cl", "Cl", 1): 242, ("Cl", "H", 1): 431,
    ("C", "F", 1): 485, ("Br", "C", 1): 276, ("Br", "Br", 1): 193,
    ("C", "S", 1): 259, ("H", "S", 1): 347, ("O", "S", 2): 523, ("S", "S", 1): 266,
    ("C", "I", 1): 240,
    # aromatic (rough, mostly cancel as spectators)
    ("C", "C", "ar"): 518, ("C", "N", "ar"): 410, ("C", "O", "ar"): 460,
}
_DEFAULT_BOND = 350.0  # fallback for unparametrised bonds


def _bond_key(b: Chem.Bond) -> tuple:
    a, c = sorted([b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()])
    bt = b.GetBondTypeAsDouble()
    order = "ar" if bt == 1.5 else int(bt)
    return (a, c, order)


def molecule_bond_energy(mol: Molecule) -> tuple[float, int]:
    """Sum of average bond enthalpies over all bonds (incl. to H).

    Returns (total_kJ, n_unknown_bonds)."""
    m = Chem.AddHs(mol.mol)
    total = 0.0
    unknown = 0
    for b in m.GetBonds():
        key = _bond_key(b)
        if key in BOND_ENTHALPY:
            total += BOND_ENTHALPY[key]
        else:
            total += _DEFAULT_BOND
            unknown += 1
    return total, unknown


def reaction_enthalpy(reactants: list[Molecule], products: list[Molecule],
                      r_coeffs: list[int], p_coeffs: list[int]) -> dict:
    e_react = unk = 0.0
    for c, m in zip(r_coeffs, reactants):
        e, u = molecule_bond_energy(m)
        e_react += c * e
        unk += u
    e_prod = 0.0
    for c, m in zip(p_coeffs, products):
        e, u = molecule_bond_energy(m)
        e_prod += c * e
        unk += u
    dH = e_react - e_prod  # kJ/mol
    return {"delta_h_kj": round(dH, 1), "unknown_bonds": int(unk)}


def thermodynamics(reactants: list[Molecule], products: list[Molecule],
                   r_coeffs: list[int] | None = None,
                   p_coeffs: list[int] | None = None,
                   temperature_k: float = 298.15) -> dict | None:
    if not reactants or not products:
        return None
    r_coeffs = r_coeffs or [1] * len(reactants)
    p_coeffs = p_coeffs or [1] * len(products)
    enth = reaction_enthalpy(reactants, products, r_coeffs, p_coeffs)
    dH = enth["delta_h_kj"]
    # crude entropy estimate from change in molecule count
    dn = sum(p_coeffs) - sum(r_coeffs)
    dS = dn * 30.0  # J/mol·K  (heuristic: more molecules → more entropy)
    dG = dH - temperature_k * dS / 1000.0  # kJ/mol
    if dG < -5:
        verdict = "خودبه‌خودی (واکنش از نظر ترمودینامیکی مطلوب است)"
    elif dG > 5:
        verdict = "غیرخودبه‌خودی (به انرژی ورودی نیاز دارد)"
    else:
        verdict = "نزدیک به تعادل (ΔG≈۰)"
    exo = "گرمازا (exothermic)" if dH < 0 else ("گرماگیر (endothermic)" if dH > 0 else "خنثی")
    return {
        "delta_h_kj": dH,
        "delta_s_j": round(dS, 1),
        "delta_g_kj": round(dG, 1),
        "temperature_k": round(temperature_k, 2),
        "type_fa": exo,
        "spontaneity_fa": verdict,
        "is_estimate": True,
        "note_fa": "تخمین به‌روش میانگین انرژی پیوند؛ آنتروپی تقریبی است.",
    }


# --- kinetics / yield ------------------------------------------------------
# base equilibrium/typical yields and whether reaction is reversible
_RULE_KINETICS = {
    "esterification": (0.66, True, "واکنش تعادلی برگشت‌پذیر؛ حذف آب بازده را بالا می‌برد"),
    "amide_formation": (0.70, True, "نیازمند گرمای بالا برای حذف آب"),
    "saponification": (0.95, False, "هیدرولیز بازی عملاً کامل پیش می‌رود"),
    "hydrogenation": (0.92, False, "افزایشی و کامل با کاتالیزور فلزی"),
    "alkene_hydration": (0.75, True, "تعادلی؛ به کاتالیزور اسیدی وابسته"),
    "alkene_halogenation": (0.95, False, "افزایش سریع و کامل"),
    "radical_halogenation": (0.5, False, "رادیکالی؛ مخلوط محصولات می‌دهد"),
    "alcohol_oxidation": (0.7, False, "به اکسنده و کنترل شرایط وابسته"),
    "neutralization": (0.99, False, "یونی و آنی؛ عملاً کامل"),
    "combustion": (0.98, False, "احتراق کامل در اکسیژن کافی"),
}


def kinetics(rule_id: str, feasible: bool, cond: Conditions) -> dict:
    base_yield, reversible, note = _RULE_KINETICS.get(rule_id, (0.7, False, "تخمین عمومی"))
    if not feasible:
        return {"feasible": False, "estimated_yield_pct": 0,
                "rate_fa": "در این شرایط انجام نمی‌شود", "relative_rate": 0.0,
                "note_fa": "شرایط لازم برای این واکنش فراهم نیست"}
    # Arrhenius-style relative rate vs 25°C (assumed Ea = 50 kJ/mol)
    Ea = 50_000.0
    T, Tref = cond.temperature_k, 298.15
    rel = math.exp(-Ea / R * (1.0 / T - 1.0 / Tref))
    # yield adjustments
    y = base_yield
    if cond.catalyst:
        y = min(0.98, y + 0.1)
    if reversible and cond.is_heated:
        y = min(0.95, y + 0.05)
    if cond.is_cold:
        y *= 0.7
    # rate class
    if rel >= 5:
        rate = "سریع"
    elif rel >= 0.5:
        rate = "متوسط"
    else:
        rate = "کند"
    if cond.catalyst:
        rate += " (کاتالیزور سرعت را افزایش می‌دهد)"
    return {
        "feasible": True,
        "estimated_yield_pct": round(y * 100),
        "rate_fa": rate,
        "relative_rate": round(rel, 2),
        "reversible": reversible,
        "note_fa": note,
        "is_estimate": True,
    }
