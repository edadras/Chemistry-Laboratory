"""
scenario.py — High-level orchestration.

A "scenario" is what the user submits: a list of reactant descriptions
(names / SMILES / formulas) plus environmental conditions, and a mode
(reaction prediction, retrosynthesis, or hypothesis/discovery). This module
parses the inputs, runs the appropriate engine, and returns a single,
structured, Persian-friendly answer including the predicted final product.
"""
from __future__ import annotations

from dataclasses import dataclass

from .conditions import Conditions
from .molecule import Molecule, MoleculeError
from .reactions import ReactionEngine
from .retrosynthesis import Retrosynthesizer
from .hypothesis import HypothesisEngine
from .pharma import DrugProfile, indication_for
from .known_compounds import resolve_name, _INDEX


def _parse_reactants(items: list[str]) -> tuple[list[Molecule], list[str]]:
    mols, errors = [], []
    for raw in items:
        raw = (raw or "").strip()
        if not raw:
            continue
        try:
            mols.append(Molecule.parse(raw))
        except MoleculeError as e:
            errors.append(f"«{raw}»: {e}")
    return mols, errors


class ScenarioProcessor:
    def __init__(self) -> None:
        self.engine = ReactionEngine()
        self.retro = Retrosynthesizer()
        self.hypo = HypothesisEngine()

    # ----- forward reaction --------------------------------------------
    def react(self, reactant_inputs: list[str], conditions: Conditions,
              pharma_mode: bool = False) -> dict:
        mols, errors = _parse_reactants(reactant_inputs)
        if not mols:
            return {"ok": False, "errors": errors or ["هیچ ماده معتبری وارد نشد"]}
        outcomes = self.engine.predict(mols, conditions, include_infeasible=True)
        feasible = [o for o in outcomes if o.feasible]
        main = feasible[0] if feasible else (outcomes[0] if outcomes else None)
        result = {
            "ok": True,
            "errors": errors,
            "reactants": [m.to_dict() for m in mols],
            "conditions": conditions.to_dict(),
            "outcomes": [o.to_dict() for o in outcomes],
            "main_product": None,
            "summary_fa": "",
        }
        if main is None:
            result["summary_fa"] = (
                "تحت این شرایط هیچ واکنش شناخته‌شده‌ای بین این مواد پیش‌بینی نشد. "
                "ممکن است مواد بی‌اثر باشند یا به شرایط/کاتالیزور دیگری نیاز باشد."
            )
            return result
        result["main_product"] = main.to_dict()
        prod_names = " + ".join(p.formula for p in main.products)
        if main.feasible:
            result["summary_fa"] = (
                f"ماده نهایی پیش‌بینی‌شده: {prod_names} از طریق «{main.rule_name_fa}». "
                f"{main.reason_fa}."
            )
        else:
            result["summary_fa"] = (
                f"محتمل‌ترین مسیر «{main.rule_name_fa}» است ولی تحت شرایط فعلی انجام نمی‌شود: "
                f"{main.reason_fa}."
            )
        if pharma_mode and main.products:
            result["pharma"] = {
                p.formula: DrugProfile(p).to_dict() for p in main.products
            }
        return result

    # ----- retrosynthesis ----------------------------------------------
    def reverse(self, target_input: str, max_depth: int = 3) -> dict:
        try:
            target = Molecule.parse(target_input)
        except MoleculeError as e:
            return {"ok": False, "errors": [str(e)]}
        steps = self.retro.analyze(target)
        blocks = self.retro.base_building_blocks(target, max_depth=max_depth)
        # pharma context if it's a known drug
        pharma = None
        name = _reverse_lookup(target.smiles)
        if name:
            ind = indication_for(name)
            pharma = {"name": name, **(ind or {})}
        summary = (
            f"هدف {target.formula} با {len(steps)} مسیر تجزیه شناسایی شد."
            if steps else
            f"برای {target.formula} مسیر تجزیه استانداردی یافت نشد؛ احتمالاً خود یک بلوک پایه است."
        )
        return {
            "ok": True,
            "target": target.to_dict(),
            "steps": [s.to_dict() for s in steps],
            "base_components": blocks,
            "pharma": pharma,
            "summary_fa": summary,
        }

    # ----- hypothesis / discovery --------------------------------------
    def hypothesize(self, reactant_inputs: list[str], conditions: Conditions,
                    pharma_mode: bool = True, top_k: int = 10) -> dict:
        mols, errors = _parse_reactants(reactant_inputs)
        if len(mols) < 2:
            return {"ok": False, "errors": (errors or []) + ["برای فرضیه‌سازی حداقل دو ماده لازم است"]}
        hyps = self.hypo.generate(mols, conditions, pharma_mode=pharma_mode, top_k=top_k)
        return {
            "ok": True,
            "errors": errors,
            "n_inputs": len(mols),
            "conditions": conditions.to_dict(),
            "hypotheses": [h.to_dict() for h in hyps],
            "summary_fa": (
                f"{len(hyps)} فرضیه تولید و آزمایش شد؛ "
                f"{sum(1 for h in hyps if h.feasible)} مورد تحت شرایط فعلی شدنی است."
            ),
        }


def _reverse_lookup(smiles: str) -> str | None:
    from rdkit import Chem
    target_canon = Chem.CanonSmiles(smiles)
    for name, smi in _INDEX.items():
        try:
            if Chem.CanonSmiles(smi) == target_canon:
                return name
        except Exception:
            continue
    return None
