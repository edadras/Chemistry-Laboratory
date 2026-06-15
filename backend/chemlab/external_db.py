"""
external_db.py — Live connectors to public chemical databases.

Clients for PubChem (PUG REST) and ChEMBL (REST). These let the system pull
real-world data — properties, bioactivities (IC50/Ki), targets — to enrich the
catalogue and to build training sets for the QSAR models.

NOTE: this execution environment restricts outbound network access by an
egress allowlist. Until `pubchem.ncbi.nlm.nih.gov` and `www.ebi.ac.uk` are
added to the environment's allowed hosts, these calls fail gracefully with a
clear message instead of raising. Once allowed, they work unchanged. ChEMBL
exports can be written straight to a CSV for `QSARModel.train_from_csv`.
"""
from __future__ import annotations

import csv
import json
import urllib.request
import urllib.parse

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def _get(url: str, timeout: int = 20) -> dict:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "chemlab/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return {"ok": True, "data": json.loads(r.read().decode())}
    except Exception as e:
        msg = str(e)
        hint = ""
        if "not in allowlist" in msg or "Forbidden" in msg or "urlopen error" in msg:
            hint = (" — احتمالاً دسترسی شبکه‌ی این محیط مسدود است؛ هاست را به "
                    "allowlist محیط اضافه کنید (تنظیمات egress).")
        return {"ok": False, "error_fa": f"خطا در اتصال: {msg}{hint}"}


# ----- PubChem ------------------------------------------------------------
def pubchem_by_name(name: str) -> dict:
    url = (f"{PUBCHEM}/compound/name/{urllib.parse.quote(name)}/property/"
           f"MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON")
    res = _get(url)
    if not res["ok"]:
        return res
    try:
        props = res["data"]["PropertyTable"]["Properties"][0]
        return {"ok": True, "source": "PubChem", "cid": props.get("CID"),
                "formula": props.get("MolecularFormula"),
                "molar_mass": props.get("MolecularWeight"),
                "smiles": props.get("CanonicalSMILES"),
                "iupac_name": props.get("IUPACName")}
    except Exception:
        return {"ok": False, "error_fa": "پاسخ PubChem قابل تفسیر نبود"}


# ----- ChEMBL -------------------------------------------------------------
def chembl_molecule(chembl_id: str) -> dict:
    res = _get(f"{CHEMBL}/molecule/{chembl_id}.json")
    if not res["ok"]:
        return res
    d = res["data"]
    structs = d.get("molecule_structures") or {}
    return {"ok": True, "source": "ChEMBL", "chembl_id": chembl_id,
            "pref_name": d.get("pref_name"),
            "smiles": structs.get("canonical_smiles"),
            "max_phase": d.get("max_phase")}


def chembl_activities(target_chembl_id: str, limit: int = 200) -> dict:
    """Fetch IC50 bioactivities for a target — ready for QSAR training."""
    url = (f"{CHEMBL}/activity.json?target_chembl_id={target_chembl_id}"
           f"&standard_type=IC50&limit={limit}")
    res = _get(url)
    if not res["ok"]:
        return res
    rows = []
    for a in res["data"].get("activities", []):
        smi = a.get("canonical_smiles")
        val = a.get("standard_value")
        if smi and val:
            rows.append({"smiles": smi, "ic50_nm": val,
                         "units": a.get("standard_units")})
    return {"ok": True, "target": target_chembl_id, "n": len(rows), "rows": rows}


def chembl_to_qsar_csv(target_chembl_id: str, out_path: str,
                       active_threshold_nm: float = 1000.0, limit: int = 500) -> dict:
    """Build a classification CSV (active if IC50 < threshold) for QSAR training."""
    act = chembl_activities(target_chembl_id, limit=limit)
    if not act["ok"]:
        return act
    written = 0
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["smiles", "label"])
        for r in act["rows"]:
            try:
                label = 1 if float(r["ic50_nm"]) < active_threshold_nm else 0
            except (TypeError, ValueError):
                continue
            w.writerow([r["smiles"], label])
            written += 1
    return {"ok": True, "path": out_path, "rows": written,
            "summary_fa": f"{written} ردیف از ChEMBL برای آموزش QSAR ذخیره شد."}
