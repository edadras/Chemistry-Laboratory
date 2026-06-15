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
from .nlp import parse_scenario
from .stoichiometry import stoichiometry


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

    # ----- natural-language scenario -----------------------------------
    def react_nl(self, text: str, pharma_mode: bool = False) -> dict:
        """Parse a free-text scenario, then run the reaction engine."""
        parsed = parse_scenario(text)
        if not parsed["reactants"]:
            return {"ok": False, "errors": ["هیچ ماده‌ی شناخته‌شده‌ای در متن پیدا نشد"],
                    "nl_parse": {"parsed_fa": parsed["parsed_fa"]}}
        res = self.react(parsed["reactants"], parsed["_conditions_obj"], pharma_mode)
        res["nl_parse"] = {"reactants": parsed["reactants"],
                           "conditions": parsed["conditions"],
                           "parsed_fa": parsed["parsed_fa"]}
        return res

    # ----- quantitative (stoichiometric) reaction ----------------------
    def react_quantitative(self, reactant_inputs: list[str], conditions: Conditions,
                           amounts: dict, actual_yield_g: float | None = None) -> dict:
        """Run a reaction and compute limiting reagent + theoretical/percent yield.

        `amounts` maps a reactant formula -> {"value": float, "unit": "g"|"mol"}.
        """
        mols, errors = _parse_reactants(reactant_inputs)
        if not mols:
            return {"ok": False, "errors": errors or ["ماده‌ی معتبری وارد نشد"]}
        outcomes = self.engine.predict(mols, conditions, include_infeasible=True)
        feasible = [o for o in outcomes if o.feasible]
        main = feasible[0] if feasible else (outcomes[0] if outcomes else None)
        if main is None:
            return {"ok": True, "main_product": None,
                    "summary_fa": "واکنشی برای محاسبه‌ی استوکیومتری یافت نشد"}
        bal = main.balanced()
        if not bal:
            return {"ok": False, "errors": ["واکنش قابل موازنه نبود؛ محاسبه‌ی کمّی ممکن نیست"]}
        amt_list = []
        for m in main.reactants_used:
            amt_list.append(amounts.get(m.formula, {"value": 1, "unit": "mol"}))
        st = stoichiometry(main.reactants_used, main.products,
                           bal["reactant_coeffs"], bal["product_coeffs"],
                           amt_list, actual_yield_g)
        return {"ok": True, "reaction": main.to_dict(include_svg=False),
                "stoichiometry": st, "summary_fa": st["summary_fa"]}

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
    def _build_route(self, target: Molecule, step, include_svg: bool = True) -> dict:
        """Turn one retro disconnection into a forward-checked synthesis route."""
        target_canon = Chem.CanonSmiles(target.smiles)
        cond, note = _SYNTH_CONDITIONS.get(step.rule_id, (Conditions(), "شرایط استاندارد"))
        confirmed = False
        balanced = thermo = kin = None
        outcomes = self.engine.predict(step.precursors, cond, include_infeasible=True)
        for oc in outcomes:
            if any(_canon_eq(p.smiles, target_canon) for p in oc.products):
                confirmed = True
                bal = oc.balanced()
                balanced = bal["equation_fa"] if bal else None
                thermo = oc.thermodynamics(bal)
                kin = oc.kinetics()
                break
        return {
            "product": target.formula,
            "product_smiles": target.smiles,
            "disconnection": step.rule_id,
            "explanation_fa": step.explanation_fa,
            "precursors": [p.to_dict(include_svg=include_svg) for p in step.precursors],
            "precursor_summary": " + ".join(p.formula for p in step.precursors),
            "conditions": cond.to_dict(),
            "conditions_note_fa": note,
            "forward_confirmed": confirmed,
            "balanced_equation": balanced,
            "thermodynamics": thermo,
            "kinetics": kin,
        }

    def synthesize(self, target_input: str, multistep: bool = False,
                   max_steps: int = 4) -> dict:
        """Plan a synthesis: target -> (retro) precursors -> (forward) confirm.

        Single-step mode lists every direct disconnection of the target.
        Multi-step mode recursively decomposes complex precursors down to
        simple, catalogue-available starting materials, returning an ordered
        recipe (base reagents → final drug), each step forward-checked.
        """
        try:
            target = Molecule.parse(target_input)
        except MoleculeError as e:
            return {"ok": False, "errors": [str(e)]}
        info = drug_info(target_input) or drug_info(_reverse_lookup(target.smiles) or "")
        base = {
            "ok": True,
            "target": target.to_dict(),
            "drug_info": info,
            "pharma": DrugProfile(target).to_dict(),
        }

        if multistep:
            steps: list[dict] = []
            leaves: list[Molecule] = []
            self._decompose(target, max_steps, steps, leaves, set())
            confirmed_n = sum(1 for s in steps if s["forward_confirmed"])
            base.update({
                "multistep": True,
                "steps": steps,
                "starting_materials": _unique_mols(leaves),
                "summary_fa": (
                    f"مسیر سنتز {len(steps)} مرحله‌ای برای {target.formula} طراحی شد؛ "
                    f"{confirmed_n} مرحله رو‌به‌جلو تأیید شد. "
                    f"مواد اولیه: {'، '.join(m['formula'] for m in _unique_mols(leaves))}."
                    if steps else
                    f"{target.formula} با قوانین فعلی تجزیه نشد؛ خود یک ماده اولیه است."
                ),
            })
            return base

        retro_steps = self.retro.analyze(target)
        routes = [self._build_route(target, s) for s in retro_steps]
        routes.sort(key=lambda r: (not r["forward_confirmed"]))
        confirmed_n = sum(1 for r in routes if r["forward_confirmed"])
        base.update({
            "multistep": False,
            "routes": routes,
            "summary_fa": (
                f"برای ساخت {target.formula}، {len(routes)} مسیر پیشنهاد شد؛ "
                f"{confirmed_n} مسیر با شبیه‌سازی رو‌به‌جلو تأیید شد (دارو بازتولید شد)."
                if routes else
                f"{target.formula} با قوانین فعلی به پیش‌سازهای ساده‌تر تجزیه نشد."
            ),
        })
        return base

    def _decompose(self, mol: Molecule, depth: int, steps: list[dict],
                   leaves: list[Molecule], visited: set[str]) -> None:
        """Post-order recursive disconnection → ordered synthesis steps."""
        from .retrosynthesis import _CLEAVAGE_RULES
        if depth <= 0 or mol.smiles in visited:
            leaves.append(mol)
            return
        visited.add(mol.smiles)
        analysis = self.retro.analyze(mol)
        cleavages = [s for s in analysis
                     if s.rule_id in _CLEAVAGE_RULES and len(s.precursors) > 1]
        if not cleavages:
            leaves.append(mol)
            return
        chosen = cleavages[0]
        # decompose complex (non-available) precursors first → deeper steps first
        for p in chosen.precursors:
            if _is_available(p):
                leaves.append(p)
            else:
                self._decompose(p, depth - 1, steps, leaves, visited)
        steps.append(self._build_route(mol, chosen, include_svg=False))

    # ----- hypothesis / discovery --------------------------------------

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


def _canon_eq(smiles_a: str, canon_b: str) -> bool:
    try:
        return Chem.CanonSmiles(smiles_a) == canon_b
    except Exception:
        return False


def _is_available(mol: Molecule) -> bool:
    """A precursor counts as a purchasable starting material if it's a known
    catalogue compound or a very small/simple molecule."""
    if _reverse_lookup(mol.smiles) is not None:
        return True
    return mol.mol.GetNumHeavyAtoms() <= 4


def _unique_mols(mols: list[Molecule]) -> list[dict]:
    seen: set[str] = set()
    out: list[dict] = []
    for m in mols:
        if m.smiles in seen:
            continue
        seen.add(m.smiles)
        name = _reverse_lookup(m.smiles)
        d = m.to_dict()
        if name:
            d["known_name"] = name
        out.append(d)
    return out


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
