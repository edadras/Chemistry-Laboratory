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

from rdkit import Chem

from .conditions import Conditions, Technique
from .molecule import Molecule, MoleculeError
from .reactions import ReactionEngine
from .retrosynthesis import Retrosynthesizer
from .hypothesis import HypothesisEngine
from .pharma import DrugProfile, indication_for
from .known_compounds import resolve_name, _INDEX, drug_info


# Suggested forward conditions per disconnection, used to turn a retro step
# into a concrete (proposed) synthesis recipe.
_SYNTH_CONDITIONS: dict[str, tuple[Conditions, str]] = {
    "ester_cut": (Conditions(temperature_c=80, catalyst="H2SO4"),
                  "استری‌شدن فیشر: کاتالیزور اسیدی (H2SO4) و گرمادهی (~۸۰°C)"),
    "amide_cut": (Conditions(temperature_c=160, techniques=[Technique.HEATING]),
                  "تراکم آمیدی: گرمادهی قوی (~۱۶۰°C) برای حذف آب"),
    "ether_cut": (Conditions(temperature_c=120, catalyst="NaOH"),
                  "سنتز ویلیامسون اتر: باز (NaOH) و گرما"),
    "reductive_amination": (Conditions(temperature_c=25, catalyst="NaBH4"),
                            "احیای آمیناسیون: تشکیل ایمین و سپس احیا با NaBH4/H2-Ni"),
    "friedel_crafts_acyl": (Conditions(temperature_c=40, catalyst="AlCl3"),
                            "آسیله‌شدن فریدل-کرافتس: کاتالیزور اسید لوویس (AlCl3)"),
    "nitro_reduction": (Conditions(temperature_c=50, catalyst="H2/Pd"),
                        "احیای گروه نیترو به آمین: H2 با کاتالیزور Pd یا Sn/HCl"),
    "alcohol_oxidation": (Conditions(temperature_c=60, catalyst="KMnO4"),
                          "اکسایش به ترکیب کربونیل با اکسنده (KMnO4/PCC)"),
}


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

    # ----- forward synthesis of a target drug --------------------------
    def synthesize(self, target_input: str) -> dict:
        """Plan a synthesis: target -> (retro) precursors -> (forward) confirm.

        For each disconnection of the target we propose precursors + suitable
        conditions, then run the forward engine to check whether the target is
        actually regenerated, yielding a verified recipe where possible.
        """
        try:
            target = Molecule.parse(target_input)
        except MoleculeError as e:
            return {"ok": False, "errors": [str(e)]}
        target_canon = Chem.CanonSmiles(target.smiles)
        info = drug_info(target_input) or drug_info(_reverse_lookup(target.smiles) or "")

        steps = self.retro.analyze(target)
        routes = []
        for s in steps:
            cond, note = _SYNTH_CONDITIONS.get(
                s.rule_id, (Conditions(), "شرایط استاندارد"))
            confirmed = False
            forward_product = None
            outcomes = self.engine.predict(s.precursors, cond, include_infeasible=True)
            for oc in outcomes:
                for p in oc.products:
                    try:
                        if Chem.CanonSmiles(p.smiles) == target_canon:
                            confirmed = True
                            forward_product = oc.to_dict(include_svg=False)
                            break
                    except Exception:
                        continue
                if confirmed:
                    break
            routes.append({
                "disconnection": s.rule_id,
                "explanation_fa": s.explanation_fa,
                "precursors": [p.to_dict() for p in s.precursors],
                "precursor_summary": " + ".join(p.formula for p in s.precursors),
                "conditions": cond.to_dict(),
                "conditions_note_fa": note,
                "forward_confirmed": confirmed,
                "forward_product": forward_product,
            })
        # confirmed routes first
        routes.sort(key=lambda r: (not r["forward_confirmed"]))
        confirmed_n = sum(1 for r in routes if r["forward_confirmed"])
        return {
            "ok": True,
            "target": target.to_dict(),
            "drug_info": info,
            "pharma": DrugProfile(target).to_dict(),
            "routes": routes,
            "summary_fa": (
                f"برای ساخت {target.formula}، {len(routes)} مسیر پیشنهاد شد؛ "
                f"{confirmed_n} مسیر با شبیه‌سازی رو‌به‌جلو تأیید شد (دارو بازتولید شد)."
                if routes else
                f"{target.formula} با قوانین فعلی به پیش‌سازهای ساده‌تر تجزیه نشد."
            ),
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
