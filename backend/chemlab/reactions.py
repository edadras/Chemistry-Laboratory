"""
reactions.py — Forward reaction engine.

Given a set of reactant molecules and `Conditions`, predicts plausible
products by:
  1. Trying each curated SMARTS reaction rule over all orderings of reactants.
  2. Applying special analytic reactions (acid-base neutralization, combustion)
     that are awkward to express as organic SMARTS.

Every predicted outcome carries the rule used, feasibility under the given
conditions, and the product molecule(s).
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from .conditions import Conditions, Technique
from .molecule import Molecule
from .reaction_data import REACTION_RULES, ReactionRule


# Pre-compile and validate all rule SMARTS once.
def _compile_rules() -> dict[str, AllChem.ChemicalReaction]:
    compiled = {}
    for rule in REACTION_RULES:
        rxn = AllChem.ReactionFromSmarts(rule.smarts)
        if rxn is None:
            raise RuntimeError(f"SMARTS نامعتبر برای قانون {rule.rid}: {rule.smarts}")
        rxn.Initialize()
        compiled[rule.rid] = rxn
    return compiled


_COMPILED = _compile_rules()


@dataclass
class Outcome:
    rule_id: str
    rule_name_fa: str
    products: list[Molecule]
    feasible: bool
    reason_fa: str
    category: str
    description_fa: str

    def to_dict(self, include_svg: bool = True) -> dict:
        return {
            "rule_id": self.rule_id,
            "rule_name_fa": self.rule_name_fa,
            "category": self.category,
            "description_fa": self.description_fa,
            "feasible": self.feasible,
            "reason_fa": self.reason_fa,
            "products": [p.to_dict(include_svg=include_svg) for p in self.products],
            "product_summary": " + ".join(p.formula for p in self.products),
        }


def _dedupe_products(prod_tuple) -> Optional[list[Molecule]]:
    mols = []
    for p in prod_tuple:
        try:
            Chem.SanitizeMol(p)
        except Exception:
            return None
        mols.append(Molecule(p))
    return mols


class ReactionEngine:
    def __init__(self) -> None:
        self.rules = {r.rid: r for r in REACTION_RULES}

    # ----- SMARTS-based prediction --------------------------------------
    def _apply_rule(self, rule: ReactionRule, reactants: list[Molecule]) -> list[list[Molecule]]:
        rxn = _COMPILED[rule.rid]
        n = rxn.GetNumReactantTemplates()
        results: list[list[Molecule]] = []
        seen_smiles: set[str] = set()
        # try every ordered selection of n reactants from the pool
        pool = [m.mol for m in reactants]
        for combo in itertools.permutations(pool, n) if len(pool) >= n else []:
            try:
                product_sets = rxn.RunReactants(combo)
            except Exception:
                continue
            for ps in product_sets:
                mols = _dedupe_products(ps)
                if mols is None:
                    continue
                key = ".".join(sorted(m.smiles for m in mols))
                if key in seen_smiles:
                    continue
                seen_smiles.add(key)
                results.append(mols)
        return results

    # ----- special inorganic reactions ----------------------------------
    def _special_reactions(self, reactants: list[Molecule], cond: Conditions) -> list[Outcome]:
        outcomes: list[Outcome] = []
        formulas = [m.formula for m in reactants]

        # Acid-base neutralization: strong acid + strong base -> salt + water
        acid = self._find(reactants, is_bronsted_acid)
        base = self._find(reactants, is_hydroxide_base)
        if acid is not None and base is not None:
            water = Molecule.from_smiles("O")
            salt = _make_salt(acid, base)
            prods = [salt, water] if salt else [water]
            outcomes.append(Outcome(
                rule_id="neutralization",
                rule_name_fa="خنثی‌سازی اسید و باز",
                products=prods,
                feasible=True,
                reason_fa="اسید و باز بلافاصله واکنش می‌دهند و نمک + آب می‌سازند",
                category="inorganic",
                description_fa="اسید + باز ⟶ نمک + آب",
            ))

        # Combustion: hydrocarbon (C,H[,O] only) + O2 -> CO2 + H2O
        o2 = self._find(reactants, lambda m: m.formula == "O2")
        fuel = self._find(reactants, is_hydrocarbon)
        if o2 is not None and fuel is not None and (cond.is_heated or cond.has(Technique.HEATING) or cond.light):
            outcomes.append(Outcome(
                rule_id="combustion",
                rule_name_fa="احتراق",
                products=[Molecule.from_smiles("O=C=O"), Molecule.from_smiles("O")],
                feasible=True,
                reason_fa="سوخت هیدروکربنی در حضور اکسیژن و جرقه/گرما می‌سوزد",
                category="inorganic",
                description_fa="هیدروکربن + O₂ ⟶ CO₂ + H₂O (+ انرژی)",
            ))
        return outcomes

    @staticmethod
    def _find(reactants, pred) -> Optional[Molecule]:
        for m in reactants:
            try:
                if pred(m):
                    return m
            except Exception:
                continue
        return None

    # ----- public API ---------------------------------------------------
    def predict(self, reactants: list[Molecule], cond: Conditions,
                include_infeasible: bool = True) -> list[Outcome]:
        outcomes: list[Outcome] = []
        for rule in REACTION_RULES:
            product_sets = self._apply_rule(rule, reactants)
            if not product_sets:
                continue
            ok, reason = rule.feasible(cond)
            if not ok and not include_infeasible:
                continue
            # take the first (most direct) product set per rule
            outcomes.append(Outcome(
                rule_id=rule.rid,
                rule_name_fa=rule.name_fa,
                products=product_sets[0],
                feasible=ok,
                reason_fa=reason,
                category=rule.category,
                description_fa=rule.description_fa,
            ))
        outcomes.extend(self._special_reactions(reactants, cond))
        # feasible first
        outcomes.sort(key=lambda o: (not o.feasible))
        return outcomes


# ----- helper chemistry predicates --------------------------------------
def is_bronsted_acid(m: Molecule) -> bool:
    # carboxylic acid OR a known strong inorganic acid
    if m.smiles in {"Cl", "OS(=O)(=O)O", "O[N+](=O)[O-]"}:
        return True
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    return m.mol.HasSubstructMatch(patt)


def is_hydroxide_base(m: Molecule) -> bool:
    return "[OH-]" in m.smiles or m.smiles in {"[Na+].[OH-]", "[K+].[OH-]", "[OH-].[Na+]"}


def is_hydrocarbon(m: Molecule) -> bool:
    syms = set(m.atom_counts().keys())
    return syms.issubset({"C", "H"}) and "C" in syms


def _make_salt(acid: Molecule, base: Molecule) -> Optional[Molecule]:
    """Best-effort salt formation for simple strong acid + strong base."""
    cation = None
    if "[Na+]" in base.smiles:
        cation = "[Na+]"
    elif "[K+]" in base.smiles:
        cation = "[K+]"
    anion_map = {
        "Cl": "[Cl-]",
        "OS(=O)(=O)O": "[O-]S(=O)(=O)[O-]",
        "O[N+](=O)[O-]": "[O-][N+](=O)[O-]",
    }
    anion = anion_map.get(acid.smiles)
    if cation and anion:
        try:
            return Molecule.from_smiles(f"{cation}.{anion}")
        except Exception:
            return None
    # carboxylic acid + NaOH -> sodium carboxylate
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if cation and acid.mol.HasSubstructMatch(patt):
        rxn = AllChem.ReactionFromSmarts("[CX3:1](=O)[OX2H1]>>[CX3:1](=O)[O-]")
        rxn.Initialize()
        try:
            prods = rxn.RunReactants((acid.mol,))
            if prods:
                p = prods[0][0]
                Chem.SanitizeMol(p)
                return Molecule.from_smiles(f"{cation}.{Chem.MolToSmiles(p)}")
        except Exception:
            return None
    return None
