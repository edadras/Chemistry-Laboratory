"""
targets.py — Biological target registry (Step 4: focus on one target).

Each entry maps a target to its ChEMBL target id and metadata so the workbench
can (a) recognise a goal like «مهارکننده COX2» and (b) pull bioactivities from
ChEMBL to train a target-specific QSAR model once network egress is allowed.
ChEMBL target ids are standard and can be edited/extended freely.
"""
from __future__ import annotations

TARGETS: dict[str, dict] = {
    "COX2": {
        "name_fa": "سیکلواکسیژناز-۲ (COX-2)",
        "chembl_id": "CHEMBL230",
        "disease_fa": "التهاب و درد (NSAIDها)",
        "aliases": ["cox2", "cox-2", "ptgs2", "سیکلواکسیژناز", "کوکس۲"],
    },
    "COX1": {
        "name_fa": "سیکلواکسیژناز-۱ (COX-1)",
        "chembl_id": "CHEMBL221",
        "disease_fa": "التهاب؛ هدف آسپرین",
        "aliases": ["cox1", "cox-1", "ptgs1"],
    },
    "EGFR": {
        "name_fa": "گیرنده‌ی فاکتور رشد اپیدرمی (EGFR)",
        "chembl_id": "CHEMBL203",
        "disease_fa": "سرطان (مهارکننده‌های تیروزین‌کیناز)",
        "aliases": ["egfr", "erbb1", "گیرنده رشد"],
    },
    "MPRO": {
        "name_fa": "پروتئاز اصلی SARS-CoV-2 (Mpro/3CLpro)",
        "chembl_id": "CHEMBL4523582",
        "disease_fa": "کووید-۱۹ (ضدویروس)",
        "aliases": ["mpro", "3clpro", "sars", "covid", "کرونا", "کووید"],
    },
    "BACE1": {
        "name_fa": "بتا-سکرتاز ۱ (BACE1)",
        "chembl_id": "CHEMBL4822",
        "disease_fa": "آلزایمر",
        "aliases": ["bace", "bace1", "آلزایمر"],
    },
}


def resolve_target(text: str | None) -> dict | None:
    """Find a target by key/alias appearing in a (possibly free-text) goal."""
    if not text:
        return None
    low = text.strip().lower()
    if low.upper() in TARGETS:
        key = low.upper()
        return {"key": key, **TARGETS[key]}
    for key, info in TARGETS.items():
        if key.lower() in low or any(a in low for a in info["aliases"]):
            return {"key": key, **info}
    return None


def list_targets() -> list[dict]:
    return [{"key": k, **v} for k, v in TARGETS.items()]
