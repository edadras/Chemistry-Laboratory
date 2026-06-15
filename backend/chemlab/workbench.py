"""
workbench.py — Drug Discovery Workbench (project v2).

A single pipeline that realises the user's vision:

    «می‌خواهم مهارکننده‌ی COX2 پیدا کنم»
        → (۱) هدف زیستی + داده‌ی ChEMBL  → (۲) مدل QSAR
        → (۳) تولید مولکول  → (۴) QSAR + ADMET + داکینگ
        → (۵) رتبه‌بندی  → (۶) مسیر سنتز (Retrosynthesis) برای هر کاندید

The differentiator (Step 7): every ranked candidate is handed to the existing
synthesis engine, so the output is not just "this molecule is good" but a
candidate *with a concrete synthetic route* — Activity + ADMET + Docking +
Retrosynthesis in one card.
"""
from __future__ import annotations

import os
import tempfile

from .qsar import QSARModel, demo_model
from .discovery import discover
from .scenario import ScenarioProcessor
from .targets import resolve_target, list_targets
from . import external_db

_PROCESSOR = ScenarioProcessor()


def _target_qsar(target: dict | None) -> tuple[QSARModel, dict]:
    """Get a QSAR model for the target: train from ChEMBL if reachable, else demo."""
    note = {"model_source": "demo", "warning_fa": None}
    if not target:
        note["warning_fa"] = "هدف مشخص نشد؛ از مدل دموی عمومی استفاده شد."
        return demo_model(), note
    # try live ChEMBL → train a target-specific model
    csv_path = os.path.join(tempfile.gettempdir(), f"qsar_{target['key']}.csv")
    built = external_db.chembl_to_qsar_csv(target["chembl_id"], csv_path, limit=400)
    if built.get("ok") and built.get("rows", 0) >= 30:
        try:
            model = QSARModel.train_from_csv(csv_path, task="classification",
                                             name=f"QSAR_{target['key']}")
            note["model_source"] = "ChEMBL"
            note["chembl_rows"] = built["rows"]
            return model, note
        except Exception as e:
            note["warning_fa"] = f"آموزش از ChEMBL ناموفق بود ({e})؛ مدل دمو استفاده شد."
    else:
        note["warning_fa"] = (
            "دسترسی به ChEMBL ممکن نشد (شبکه‌ی محیط مسدود است)؛ از مدل دموی "
            "BBB استفاده شد. پس از افزودن www.ebi.ac.uk به allowlist، مدل "
            "مخصوص همین هدف از داده‌ی واقعی آموزش می‌بیند.")
    return demo_model(), note


def run(goal: str | None = None, target_key: str | None = None,
        population_size: int = 40, generations: int = 8, top_k: int = 8,
        use_docking: bool = False, with_routes: bool = True) -> dict:
    """Run the full workbench pipeline and return ranked candidates with routes."""
    target = resolve_target(target_key) or resolve_target(goal)
    qsar, note = _target_qsar(target)

    disc = discover(qsar=qsar, objective="activity",
                    population_size=population_size, generations=generations,
                    top_k=top_k, use_docking=use_docking)

    # Step 7 — attach a synthetic route to every candidate.
    if with_routes:
        for c in disc.get("candidates", []):
            plan = _PROCESSOR.synthesize(c["smiles"], multistep=True, max_steps=4)
            steps = plan.get("steps", []) if plan.get("ok") else []
            c["retrosynthesis"] = {
                "n_steps": len(steps),
                "steps": [{"product": s["product"],
                           "precursors": s["precursor_summary"],
                           "disconnection": s["disconnection"],
                           "conditions_note_fa": s["conditions_note_fa"],
                           "forward_confirmed": s["forward_confirmed"]}
                          for s in steps],
                "starting_materials": [m["formula"] for m in plan.get("starting_materials", [])],
                "summary_fa": (plan.get("summary_fa", "")
                               if steps else
                               "مسیر سنتز استانداردی یافت نشد؛ این مولکول نو احتمالاً "
                               "به طراحی سنتز اختصاصی نیاز دارد."),
            }

    stages = [
        "۱) شناسایی هدف زیستی" + (f": {target['name_fa']}" if target else " (عمومی)"),
        f"۲) مدل QSAR ({note['model_source']})",
        f"۳) تولید مولکول‌های نو ({disc.get('molecules_evaluated', 0)} مولکول)",
        "۴) ارزیابی QSAR + ADMET" + ("+ داکینگ" if use_docking else ""),
        "۵) رتبه‌بندی با امتیاز ترکیبی",
        "۶) طراحی مسیر سنتز برای هر کاندید" if with_routes else "۶) (مسیر سنتز غیرفعال)",
    ]
    return {
        "ok": True,
        "goal": goal,
        "target": target,
        "qsar_model": qsar.to_dict(),
        "model_note": note,
        "pipeline_stages_fa": stages,
        "candidates": disc.get("candidates", []),
        "molecules_evaluated": disc.get("molecules_evaluated", 0),
        "summary_fa": (
            f"میز کار برای هدف «{target['name_fa'] if target else 'عمومی'}» اجرا شد؛ "
            f"{disc.get('molecules_evaluated', 0)} مولکول تولید و {len(disc.get('candidates', []))} "
            f"کاندیدای برتر با فعالیت/ADMET" + ("/داکینگ" if use_docking else "") +
            " و مسیر سنتز ارائه شد."
        ),
        "available_targets": list_targets(),
    }
