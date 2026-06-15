"""
pharma.py — Pharmaceutical / drug-discovery analysis.

Focused on the treatment & medicine domain the user cares about. Provides:
  * Drug-likeness assessment (Lipinski's Rule of Five, Veber rules).
  * A crude QED-style desirability score.
  * Structural-alert screening for a few classic toxicophores.
  * A small knowledge base linking known drugs to indications (Persian).
"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski, rdMolDescriptors

from .molecule import Molecule
from .known_compounds import drug_info

# (SMARTS, name_fa) toxicophore alerts — illustrative, not exhaustive.
STRUCTURAL_ALERTS = [
    ("[N+](=O)[O-]", "گروه نیترو (سمیت بالقوه)"),
    ("N=[N+]=[N-]", "آزید (ناپایدار/سمی)"),
    ("[CX3](=O)[OX2][CX3](=O)", "انیدرید (واکنش‌پذیر)"),
    ("[#6][Cl,Br,I]", "هالید آلکیلی (آلکیله‌کننده بالقوه)"),
    ("[SX2H]", "تیول آزاد"),
]


@dataclass
class DrugProfile:
    molecule: Molecule

    def lipinski(self) -> dict:
        m = self.molecule.mol
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        hbd = Lipinski.NumHDonors(m)
        hba = Lipinski.NumHAcceptors(m)
        violations = []
        if mw > 500:
            violations.append("جرم مولکولی > ۵۰۰")
        if logp > 5:
            violations.append("LogP > ۵ (چربی‌دوستی زیاد)")
        if hbd > 5:
            violations.append("دهنده پیوند هیدروژنی > ۵")
        if hba > 10:
            violations.append("گیرنده پیوند هیدروژنی > ۱۰")
        return {
            "molecular_weight": round(mw, 2),
            "logp": round(logp, 2),
            "h_bond_donors": hbd,
            "h_bond_acceptors": hba,
            "violations": violations,
            "passes": len(violations) <= 1,  # RO5: at most one violation
        }

    def veber(self) -> dict:
        m = self.molecule.mol
        rotb = Lipinski.NumRotatableBonds(m)
        tpsa = rdMolDescriptors.CalcTPSA(m)
        ok = rotb <= 10 and tpsa <= 140
        return {"rotatable_bonds": rotb, "tpsa": round(tpsa, 1),
                "oral_bioavailability_ok": ok}

    def qed(self) -> float:
        try:
            return round(QED.qed(self.molecule.mol), 3)
        except Exception:
            return float("nan")

    def alerts(self) -> list[str]:
        hits = []
        for smarts, label in STRUCTURAL_ALERTS:
            patt = Chem.MolFromSmarts(smarts)
            if patt is not None and self.molecule.mol.HasSubstructMatch(patt):
                hits.append(label)
        return hits

    def verdict_fa(self) -> str:
        ro5 = self.lipinski()
        veber = self.veber()
        qed = self.qed()
        if ro5["passes"] and veber["oral_bioavailability_ok"] and qed >= 0.5:
            return "پروفایل دارویی مطلوب: کاندیدای امیدوارکننده برای داروی خوراکی"
        if ro5["passes"]:
            return "از قانون لیپینسکی عبور می‌کند ولی برخی شاخص‌ها بهینه نیستند"
        return "نقض قوانین دارورسانی خوراکی — به‌احتمال جذب خوراکی ضعیف"

    def to_dict(self) -> dict:
        return {
            "lipinski": self.lipinski(),
            "veber": self.veber(),
            "qed": self.qed(),
            "structural_alerts": self.alerts(),
            "verdict_fa": self.verdict_fa(),
        }


def indication_for(name: str) -> dict | None:
    info = drug_info(name)
    if not info:
        return None
    return {"indication_fa": info.get("indication_fa"),
            "drug_class_fa": info.get("drug_class_fa"),
            "category": info.get("category")}
