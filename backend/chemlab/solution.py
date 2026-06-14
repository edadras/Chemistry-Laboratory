"""
solution.py — Solution chemistry: pH, pKa, acid–base equilibria, solubility.

Computes the pH of an aqueous solution of an acid or base at a given molar
concentration using standard weak/strong equilibrium approximations, looks up
pKa for common species, and gives a heuristic aqueous-solubility verdict from
logP / polarity descriptors.
"""
from __future__ import annotations

import math

from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors

from .molecule import Molecule

Kw = 1.0e-14

# pKa of common acids (first dissociation) and conjugate-acid pKa of bases.
# Keyed by canonical SMILES.
PKA: dict[str, dict] = {
    "Cl":             {"pka": -7.0, "type": "acid", "name_fa": "اسید کلریدریک (قوی)"},
    "OS(=O)(=O)O":    {"pka": -3.0, "type": "acid", "name_fa": "اسید سولفوریک (قوی)"},
    "O=[N+]([O-])O":  {"pka": -1.4, "type": "acid", "name_fa": "اسید نیتریک (قوی)"},
    "CC(=O)O":        {"pka": 4.76, "type": "acid", "name_fa": "اسید استیک (ضعیف)"},
    "O=CO":           {"pka": 3.75, "type": "acid", "name_fa": "اسید فرمیک (ضعیف)"},
    "O=C(O)c1ccccc1O":{"pka": 2.97, "type": "acid", "name_fa": "سالیسیلیک اسید"},
    "CC(=O)Oc1ccccc1C(=O)O": {"pka": 3.49, "type": "acid", "name_fa": "آسپرین"},
    "Oc1ccccc1":      {"pka": 9.95, "type": "acid", "name_fa": "فنل (بسیار ضعیف)"},
    "[OH-]":          {"pka": 15.7, "type": "strong_base", "name_fa": "هیدروکسید"},
    "N":              {"pka": 9.25, "type": "base", "name_fa": "آمونیاک (پایه‌ی ضعیف)"},
    "C":              {"pka": 50.0, "type": "neutral", "name_fa": "متان (خنثی)"},
}

STRONG_ACID_SMILES = {"Cl", "OS(=O)(=O)O", "O=[N+]([O-])O"}
STRONG_BASE_HINTS = ("[OH-]", "[Na+].[OH-]", "[K+].[OH-]")


def lookup_pka(mol: Molecule) -> dict | None:
    try:
        canon = Chem.CanonSmiles(mol.smiles)
    except Exception:
        return None
    return PKA.get(canon)


def ph_of_solution(mol: Molecule, concentration_m: float) -> dict:
    """Estimate the pH of an aqueous solution at molar concentration C."""
    if concentration_m <= 0:
        return {"ok": False, "error_fa": "غلظت باید مثبت باشد"}
    info = lookup_pka(mol)
    canon = Chem.CanonSmiles(mol.smiles)
    C = concentration_m

    # strong base (hydroxide salts)
    if any(h in mol.smiles for h in STRONG_BASE_HINTS):
        oh = C
        pOH = -math.log10(oh)
        ph = 14 - pOH
        kind = "باز قوی"
    elif canon in STRONG_ACID_SMILES:
        ph = -math.log10(C)
        kind = "اسید قوی"
    elif info and info["type"] == "acid":
        Ka = 10 ** (-info["pka"])
        # weak acid: [H+] = sqrt(Ka*C) (approx)
        h = math.sqrt(Ka * C)
        ph = -math.log10(h)
        kind = "اسید ضعیف"
    elif info and info["type"] == "base":
        Kb = Kw / 10 ** (-info["pka"])
        oh = math.sqrt(Kb * C)
        ph = 14 + math.log10(oh)
        kind = "پایه‌ی ضعیف"
    else:
        ph = 7.0
        kind = "خنثی (یا pKa نامشخص)"

    ph = max(0.0, min(14.0, ph))
    if ph < 6.5:
        nature = "اسیدی"
    elif ph > 7.5:
        nature = "بازی"
    else:
        nature = "خنثی"
    return {
        "ok": True,
        "ph": round(ph, 2),
        "concentration_m": C,
        "kind_fa": kind,
        "nature_fa": nature,
        "pka": info["pka"] if info else None,
        "name_fa": info["name_fa"] if info else None,
        "summary_fa": f"محلول {C} مولار {info['name_fa'] if info else mol.formula}: "
                      f"pH ≈ {round(ph,2)} ({nature}، {kind}).",
    }


def solubility(mol: Molecule) -> dict:
    """Heuristic aqueous solubility from logP and polar surface area."""
    logp = Crippen.MolLogP(mol.mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol.mol)
    mw = rdMolDescriptors.CalcExactMolWt(mol.mol)
    # crude ESOL-like trend: more lipophilic & heavier → less soluble
    if logp < 0 and tpsa > 40:
        verdict = "بسیار محلول در آب (قطبی)"
    elif logp < 1.5:
        verdict = "محلول در آب"
    elif logp < 3:
        verdict = "کم‌محلول در آب، محلول در حلال آلی"
    else:
        verdict = "تقریباً نامحلول در آب (چربی‌دوست)"
    return {"logp": round(logp, 2), "tpsa": round(tpsa, 1),
            "molar_mass": round(mw, 2), "verdict_fa": verdict,
            "like_dissolves_like_fa": "آب‌دوست" if logp < 1 else "چربی‌دوست"}
