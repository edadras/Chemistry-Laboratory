"""
docking.py — Molecular docking with AutoDock Vina (Phase 4).

Pipeline: SMILES → 3D embed (RDKit ETKDG + MMFF) → PDBQT (Meeko) → Vina dock
against a rigid receptor → binding affinity (kcal/mol, more negative = stronger).

A small *demo receptor* (a peptide scaffold) is generated and cached so the
whole pipeline runs end-to-end out of the box. It is only a smoke-test target,
NOT a real protein — for real work, point `receptor_pdbqt` at a prepared
receptor (e.g. from the PDB via Meeko's `mk_prepare_receptor`). If Vina/Meeko
are unavailable, a transparent physicochemical estimate is returned instead.
"""
from __future__ import annotations

import os
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors

_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
_DEMO_RECEPTOR = os.path.join(_DATA_DIR, "demo_receptor.pdbqt")

try:
    from vina import Vina
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    _HAS_VINA = True
except Exception:
    _HAS_VINA = False


def _embed3d(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol


def _to_pdbqt(mol: Chem.Mol) -> str | None:
    try:
        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        out = PDBQTWriterLegacy.write_string(setups[0])
        # meeko returns (string, is_ok, err) in recent versions
        return out[0] if isinstance(out, tuple) else out
    except Exception:
        return None


def _rigidify(pdbqt: str) -> str:
    """Keep only atom records → a valid rigid-receptor PDBQT (no torsion tree)."""
    keep = [ln for ln in pdbqt.splitlines()
            if ln.startswith(("ATOM", "HETATM"))]
    return "\n".join(keep) + "\n"


def _coords_centroid(pdbqt: str) -> tuple[float, float, float]:
    xs, ys, zs = [], [], []
    for ln in pdbqt.splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            xs.append(float(ln[30:38])); ys.append(float(ln[38:46])); zs.append(float(ln[46:54]))
    n = max(1, len(xs))
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def ensure_demo_receptor() -> str | None:
    """Build & cache a small demo receptor PDBQT (peptide scaffold)."""
    if os.path.exists(_DEMO_RECEPTOR):
        return _DEMO_RECEPTOR
    if not _HAS_VINA:
        return None
    # a short peptide-like molecule as a stand-in binding pocket
    pep = "CC(C)C[C@@H](C(=O)O)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CO"
    mol = _embed3d(pep)
    if mol is None:
        return None
    pdbqt = _to_pdbqt(mol)
    if pdbqt is None:
        return None
    os.makedirs(_DATA_DIR, exist_ok=True)
    with open(_DEMO_RECEPTOR, "w") as f:
        f.write(_rigidify(pdbqt))
    return _DEMO_RECEPTOR


def _estimate_affinity(smiles: str) -> dict:
    """Fallback when Vina is unavailable: crude affinity proxy from logP/MW."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"ok": False, "error_fa": "SMILES نامعتبر"}
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    aff = -(2.0 + 0.6 * logp + 0.004 * mw)  # heuristic kcal/mol
    return {"ok": True, "engine": "estimate",
            "binding_affinity_kcal_mol": round(aff, 2),
            "is_estimate": True,
            "note_fa": "تخمین فیزیکوشیمیایی (Vina در دسترس نبود)."}


def dock(smiles: str, receptor_pdbqt: str | None = None,
         center: tuple[float, float, float] | None = None,
         box_size: tuple[float, float, float] = (22, 22, 22),
         exhaustiveness: int = 4) -> dict:
    """Dock a ligand and return its best binding affinity (kcal/mol)."""
    if not _HAS_VINA:
        return _estimate_affinity(smiles)
    receptor = receptor_pdbqt or ensure_demo_receptor()
    if receptor is None or not os.path.exists(receptor):
        return _estimate_affinity(smiles)

    mol = _embed3d(smiles)
    if mol is None:
        return {"ok": False, "error_fa": "ساخت ساختار سه‌بعدی لیگاند ناموفق بود"}
    lig_pdbqt = _to_pdbqt(mol)
    if lig_pdbqt is None:
        return {"ok": False, "error_fa": "تبدیل لیگاند به PDBQT ناموفق بود"}

    if center is None:
        with open(receptor) as f:
            center = _coords_centroid(f.read())
    try:
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(rigid_pdbqt_filename=receptor)
        v.set_ligand_from_string(lig_pdbqt)
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=5)
        energies = v.energies(n_poses=1)
        best = float(energies[0][0])
        return {"ok": True, "engine": "AutoDock Vina",
                "binding_affinity_kcal_mol": round(best, 2),
                "receptor": os.path.basename(receptor),
                "demo_receptor": receptor == _DEMO_RECEPTOR,
                "is_estimate": False,
                "note_fa": ("اتصال به گیرنده‌ی نمونه (آزمایشی) — برای کار واقعی "
                            "گیرنده‌ی پروتئینی آماده بدهید."
                            if receptor == _DEMO_RECEPTOR else
                            "داکینگ با گیرنده‌ی ارائه‌شده.")}
    except Exception as e:
        return {"ok": False, "error_fa": f"خطا در اجرای Vina: {str(e)[:120]}",
                "fallback": _estimate_affinity(smiles)}
