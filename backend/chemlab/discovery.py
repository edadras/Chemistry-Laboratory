"""
discovery.py — Closed-loop de novo drug discovery & candidate ranking.

Ties the pieces together (the user's Phase 7):
    generate (GA) → predict activity (QSAR) → ADMET → synthesizability → rank

Given a trained QSAR model and an objective, it evolves novel molecules and
returns a ranked candidate list, each with Activity / ADMET / synthesizability
sub-scores and a composite score — exactly the "Candidate #1 / #2" output.
"""
from __future__ import annotations

import os
import sys

from rdkit import Chem
from rdkit.Chem import QED, RDConfig

from .molecule import Molecule
from .qsar import QSARModel, demo_model
from .generator import evolve
from .admet import admet

# Synthetic-accessibility score (RDKit Contrib sascorer, 1=easy … 10=hard).
try:
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer  # type: ignore
    _HAS_SA = True
except Exception:
    _HAS_SA = False


def synth_accessibility(mol: Chem.Mol) -> float:
    if _HAS_SA:
        try:
            return float(sascorer.calculateScore(mol))
        except Exception:
            pass
    # fallback heuristic: more rings/stereo/atoms → harder
    nrings = mol.GetRingInfo().NumRings()
    natoms = mol.GetNumHeavyAtoms()
    nstereo = len(Chem.FindMolChiralCenters(mol, useLegacyImplementation=False))
    return min(10.0, 1.0 + 0.05 * natoms + 0.5 * nrings + 0.4 * nstereo)


def _sa_to_score(sa: float) -> float:
    """Map SA (1 easy … 10 hard) → 0..1 (1 = easy to synthesise)."""
    return max(0.0, min(1.0, (10.0 - sa) / 9.0))


def candidate_profile(smiles: str, qsar: QSARModel) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    pred = qsar.predict(smiles)
    activity = pred.get("score", 0.0)
    if qsar.task == "regression":  # squash regression into 0..1 for ranking
        activity = max(0.0, min(1.0, activity / 10.0))
    qed = float(QED.qed(mol))
    sa = synth_accessibility(mol)
    sa_score = _sa_to_score(sa)
    adm = admet(Molecule(mol))
    admet_score = _admet_score(adm)
    composite = round(100 * (0.45 * activity + 0.30 * admet_score +
                             0.25 * sa_score), 1)
    return {
        "smiles": Chem.MolToSmiles(mol),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "activity_score": round(100 * activity, 1),
        "admet_score": round(100 * admet_score, 1),
        "synthesizability_score": round(100 * sa_score, 1),
        "sa_raw": round(sa, 2),
        "qed": round(qed, 3),
        "composite_score": composite,
        "admet": adm,
        "qsar_prediction": pred,
    }


def _admet_score(adm: dict) -> float:
    s = 0.0
    if adm.get("gi_absorption_high"):
        s += 0.5
    if not adm.get("toxicophores"):
        s += 0.3
    d = adm.get("descriptors", {})
    if d.get("mw", 999) <= 500 and d.get("logp", 9) <= 5:
        s += 0.2
    return min(1.0, s)


def discover(qsar: QSARModel | None = None, objective: str = "activity",
             seeds: list[str] | None = None, population_size: int = 40,
             generations: int = 8, top_k: int = 10) -> dict:
    """Run the full generate→score→rank loop and return ranked candidates."""
    qsar = qsar or demo_model()

    def fitness(smi: str) -> float:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return -1.0
        pred = qsar.predict(smi)
        act = pred.get("score", 0.0)
        if qsar.task == "regression":
            act = max(0.0, min(1.0, act / 10.0))
        qed = float(QED.qed(mol))
        sa_score = _sa_to_score(synth_accessibility(mol))
        if objective == "qed":
            return 0.7 * qed + 0.3 * sa_score
        # default: activity-driven, regularised by drug-likeness & synthesizability
        return 0.55 * act + 0.25 * qed + 0.20 * sa_score

    result = evolve(fitness, seeds=seeds, population_size=population_size,
                    generations=generations)
    candidates = []
    for smi, _fit in result.population[:top_k]:
        prof = candidate_profile(smi, qsar)
        if prof:
            candidates.append(prof)
    candidates.sort(key=lambda c: c["composite_score"], reverse=True)
    for i, c in enumerate(candidates, 1):
        c["rank"] = i
    return {
        "ok": True,
        "model": qsar.to_dict(),
        "objective": objective,
        "generations": result.generations,
        "molecules_evaluated": result.n_evaluated,
        "candidates": candidates,
        "summary_fa": (
            f"{result.n_evaluated} مولکول نو تولید و ارزیابی شد؛ "
            f"{len(candidates)} کاندیدای برتر بر اساس امتیاز ترکیبی "
            f"(فعالیت + ADMET + سنتزپذیری) رتبه‌بندی شدند."
        ),
    }
