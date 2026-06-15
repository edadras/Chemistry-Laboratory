"""
stoichiometry.py — Quantitative stoichiometry.

Given a balanced reaction and the amounts of reactants (in grams or moles),
computes moles of each species, identifies the limiting reagent, and reports
theoretical yield of each product. If an actual yield is supplied, percent
yield is computed too.
"""
from __future__ import annotations

from .molecule import Molecule


def _to_moles(value: float, unit: str, molar_mass: float) -> float:
    unit = (unit or "g").lower()
    if unit in ("mol", "mole", "moles"):
        return value
    if unit in ("g", "gram", "grams"):
        return value / molar_mass
    if unit in ("mg",):
        return value / 1000.0 / molar_mass
    if unit in ("kg",):
        return value * 1000.0 / molar_mass
    raise ValueError(f"واحد ناشناخته: {unit}")


def stoichiometry(reactants: list[Molecule], products: list[Molecule],
                  r_coeffs: list[int], p_coeffs: list[int],
                  amounts: list[dict], actual_yield_g: float | None = None,
                  product_index: int = 0) -> dict:
    """Quantitative analysis of a balanced reaction.

    amounts: list aligned with reactants, each {"value": float, "unit": "g"|"mol"}.
    """
    r_moles = []
    for m, amt in zip(reactants, amounts):
        val = float(amt.get("value", 0))
        moles = _to_moles(val, amt.get("unit", "g"), m.molar_mass)
        r_moles.append(moles)

    # limiting reagent = min(moles_i / coeff_i)
    ratios = [(r_moles[i] / r_coeffs[i] if r_coeffs[i] else float("inf"))
              for i in range(len(reactants))]
    limiting_idx = min(range(len(ratios)), key=lambda i: ratios[i])
    limiting_ratio = ratios[limiting_idx]

    reactant_rows = []
    for i, m in enumerate(reactants):
        consumed = limiting_ratio * r_coeffs[i]
        reactant_rows.append({
            "formula": m.formula,
            "coeff": r_coeffs[i],
            "supplied_mol": round(r_moles[i], 4),
            "supplied_g": round(r_moles[i] * m.molar_mass, 3),
            "consumed_mol": round(consumed, 4),
            "leftover_mol": round(max(0.0, r_moles[i] - consumed), 4),
            "is_limiting": i == limiting_idx,
        })

    product_rows = []
    for j, m in enumerate(products):
        moles = limiting_ratio * p_coeffs[j]
        product_rows.append({
            "formula": m.formula,
            "coeff": p_coeffs[j],
            "theoretical_mol": round(moles, 4),
            "theoretical_g": round(moles * m.molar_mass, 3),
        })

    result = {
        "ok": True,
        "limiting_reagent": reactants[limiting_idx].formula,
        "reactants": reactant_rows,
        "products": product_rows,
        "summary_fa": (
            f"واکنش‌دهنده‌ی محدودکننده: {reactants[limiting_idx].formula}. "
            f"بازده تئوری {products[product_index].formula}: "
            f"{product_rows[product_index]['theoretical_g']} گرم."
        ),
    }
    if actual_yield_g is not None and product_rows:
        theo = product_rows[product_index]["theoretical_g"]
        pct = (actual_yield_g / theo * 100) if theo else 0
        result["percent_yield"] = round(pct, 1)
        result["summary_fa"] += f" بازده واقعی: {round(pct,1)}%."
    return result
