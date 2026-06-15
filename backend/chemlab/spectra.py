"""
spectra.py — Predicted spectroscopy: IR, ¹H NMR and mass spectrometry.

Educational, rule-based predictions:
  * IR — characteristic absorption bands from functional-group SMARTS.
  * ¹H NMR — distinct proton environments (RDKit symmetry classes) with rough
    chemical-shift ranges and integration (H count).
  * MS — molecular-ion m/z (and [M+H]⁺) plus a few common neutral losses.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .molecule import Molecule

# (SMARTS, wavenumber band cm⁻¹, assignment_fa)
IR_BANDS = [
    ("[OX2H][CX3]=O", "2500–3300 (پهن)", "O–H کربوکسیلیک اسید"),
    ("[OX2H]", "3200–3550 (پهن)", "کشش O–H (الکل/فنل)"),
    ("[NX3;H1,H2]", "3300–3500", "کشش N–H (آمین/آمید)"),
    ("[CX3](=O)[OX2H]", "1700–1725", "C=O کربوکسیلیک اسید"),
    ("[CX3](=O)[OX2][#6]", "1735–1750", "C=O استر"),
    ("[CX3](=O)[NX3]", "1630–1700", "C=O آمید"),
    ("[CX3H1](=O)", "1720–1740", "C=O آلدهید"),
    ("[#6][CX3](=O)[#6]", "1705–1725", "C=O کتون"),
    ("[CX2]#[CX2]", "2100–2260", "کشش C≡C"),
    ("[CX2]#[NX1]", "2210–2260", "کشش C≡N (نیتریل)"),
    ("[$([NX3](=O)=O),$([NX3+](=O)[O-])]", "1500 و 1350", "کشش نامتقارن/متقارن NO₂"),
    ("c1ccccc1", "1450–1600", "کشش C=C آروماتیک"),
    ("[CX3]=[CX3]", "1620–1680", "کشش C=C آلکن"),
    ("[CX4H]", "2850–2960", "کشش C–H اشباع"),
]


def ir_spectrum(mol: Molecule) -> dict:
    bands = []
    for smarts, band, assign in IR_BANDS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.mol.HasSubstructMatch(patt):
            bands.append({"wavenumber_cm": band, "assignment_fa": assign})
    return {"bands": bands,
            "summary_fa": (f"{len(bands)} نوار جذبی شاخص پیش‌بینی شد."
                           if bands else "نوار شاخصی شناسایی نشد.")}


def nmr_h_spectrum(mol: Molecule) -> dict:
    m = Chem.AddHs(mol.mol)
    # symmetry classes group equivalent H atoms into one signal
    ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    env_groups: dict[int, dict] = {}
    for atom in m.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        nbrs = atom.GetNeighbors()
        heavy_idx = nbrs[0].GetIdx() if nbrs else -1
        cls = ranks[atom.GetIdx()]
        env_groups.setdefault(cls, {"count": 0, "heavy_idx": heavy_idx})
        env_groups[cls]["count"] += 1

    signals = []
    for info in env_groups.values():
        shift, label, mid = _classify_h(m, info["heavy_idx"])
        signals.append({"shift_ppm": shift, "integration_h": info["count"],
                        "assignment_fa": label, "_mid": mid})
    signals.sort(key=lambda s: s["_mid"])
    for s in signals:
        s.pop("_mid", None)
    return {"num_signals": len(signals), "signals": signals,
            "summary_fa": f"{len(signals)} محیط پروتونی غیرمعادل پیش‌بینی شد."}


def _classify_h(m: Chem.Mol, heavy_idx: int):
    """Return (shift_range_fa, assignment_fa, sort_midpoint) for an H on heavy_idx."""
    if heavy_idx < 0:
        return ("1–5", "H متصل به هترواتم", 3.0)
    a = m.GetAtomWithIdx(heavy_idx)
    sym = a.GetSymbol()
    if sym == "O":
        # carboxylic acid OH vs alcohol/phenol OH
        for nb in a.GetNeighbors():
            if nb.GetSymbol() == "C" and any(
                    b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(nb).GetSymbol() == "O"
                    for b in nb.GetBonds()):
                return ("10–12", "H کربوکسیلیک اسید", 11.0)
        return ("1–5 (متغیر)", "H هیدروکسیل/فنل", 4.0)
    if sym == "N":
        return ("1–8 (متغیر)", "H آمین/آمید", 4.5)
    if sym == "C":
        if a.GetIsAromatic():
            return ("6.5–8.5", "H آروماتیک", 7.5)
        # aldehyde: C(=O)H
        if any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(a).GetSymbol() == "O"
               for b in a.GetBonds()):
            return ("9–10", "H آلدهید", 9.5)
        nbr_syms = [n.GetSymbol() for n in a.GetNeighbors()]
        # alpha to a carbonyl
        for n in a.GetNeighbors():
            if n.GetSymbol() == "C" and any(
                    b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(n).GetSymbol() == "O"
                    for b in n.GetBonds()):
                return ("2.0–2.6", "H آلفای کربونیل", 2.3)
        if "O" in nbr_syms:
            return ("3.3–4.5", "H روی کربن متصل به اکسیژن", 3.9)
        if "N" in nbr_syms:
            return ("2.2–3.5", "H روی کربن متصل به نیتروژن", 2.8)
        if a.GetHybridization() == Chem.HybridizationType.SP2:
            return ("4.5–6.5", "H وینیلی (آلکن)", 5.5)
        return ("0.8–1.8", "H آلکیل اشباع", 1.3)
    return ("1–5", "H نامشخص", 3.0)


def mass_spectrum(mol: Molecule) -> dict:
    exact = Descriptors.ExactMolWt(mol.mol)
    avg = Descriptors.MolWt(mol.mol)
    losses = [("[M]⁺", 0, "یون مولکولی"),
              ("[M+H]⁺", 1.0078, "پروتونه‌شده"),
              ("[M−H₂O]⁺", -18.0106, "حذف آب"),
              ("[M−CO]⁺", -27.9949, "حذف CO"),
              ("[M−CH₃]⁺", -15.0235, "حذف متیل"),
              ("[M−OH]⁺", -17.0027, "حذف هیدروکسیل")]
    peaks = [{"label": lbl, "mz": round(exact + d, 4), "note_fa": note}
             for lbl, d, note in losses]
    return {"exact_mass": round(exact, 4), "average_mass": round(avg, 3),
            "molecular_ion_mz": round(exact, 4), "peaks": peaks,
            "summary_fa": f"یون مولکولی [M]⁺ در m/z ≈ {round(exact,2)}."}


def full_spectra(mol: Molecule) -> dict:
    return {"ir": ir_spectrum(mol), "nmr_h": nmr_h_spectrum(mol),
            "ms": mass_spectrum(mol)}
