"""
retrosynthesis.py — Reverse engineering of a target molecule.

Given a final product (formula → molecule), work backwards to plausible
precursors. Two complementary strategies are used:

  1. Rule reversal: every forward bond-forming rule is reversed
     (`products >> reactants`) and applied to the target, disconnecting it
     into simpler building blocks.
  2. Functional-group disconnection: generic disconnections (ester, amide)
     that map a target onto classic starting materials.

The result is a (shallow) retro tree: each node is a set of precursors that
would react to form the target under stated conditions.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import AllChem

from .molecule import Molecule
from .reaction_data import REACTION_RULES


def _reverse_smarts(smarts: str) -> str:
    left, right = smarts.split(">>")
    return f"{right}>>{left}"


# Reactions that form a new skeleton (worth reversing). Pure additions of small
# molecules (hydration etc.) are reversible too but produce trivial precursors.
_RETRO_RULES = {
    "esterification": "استر را به کربوکسیلیک‌اسید + الکل تجزیه می‌کند",
    "amide_formation": "آمید را به کربوکسیلیک‌اسید + آمین تجزیه می‌کند",
}

# Dedicated, more permissive retro disconnections (independent of forward rules).
# Using [#6] matches both aromatic and aliphatic carbons, so aryl esters
# (e.g. aspirin's acetate) are disconnected correctly.
_GENERIC_DISCONNECTIONS = [
    ("ester_cut", "برش پیوند استری: استر ⟶ کربوکسیلیک‌اسید + الکل/فنل",
     "[C:1](=[O:2])[O:3][#6:4]>>[C:1](=[O:2])[OX2H:3].[#6:4][OX2H]"),
    ("amide_cut", "برش پیوند آمیدی: آمید ⟶ کربوکسیلیک‌اسید + آمین",
     "[C:1](=[O:2])[NX3:3]>>[C:1](=[O:2])[OX2H].[NX3:3][H]"),
    ("ether_cut", "برش پیوند اتری: اتر ⟶ دو الکل/فنل",
     "[#6:1][OX2;!$(O[CX3]=[OX1]):2][#6:3]>>[#6:1][OX2:2][H].[#6:3][OX2H]"),
]


@dataclass
class RetroStep:
    rule_id: str
    explanation_fa: str
    precursors: list[Molecule]

    def to_dict(self, include_svg: bool = True) -> dict:
        return {
            "rule_id": self.rule_id,
            "explanation_fa": self.explanation_fa,
            "precursors": [p.to_dict(include_svg=include_svg) for p in self.precursors],
            "precursor_summary": " + ".join(p.formula for p in self.precursors),
        }


class Retrosynthesizer:
    def __init__(self) -> None:
        self._rules = {r.rid: r for r in REACTION_RULES}
        self._reverse = {}
        for rid, expl in _RETRO_RULES.items():
            rule = self._rules[rid]
            rxn = AllChem.ReactionFromSmarts(_reverse_smarts(rule.smarts))
            if rxn is not None:
                rxn.Initialize()
                self._reverse[rid] = (rxn, expl, rule.name_fa)
        # generic disconnections (expl is "label: detail")
        for rid, expl, smarts in _GENERIC_DISCONNECTIONS:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is not None:
                rxn.Initialize()
                label, _, detail = expl.partition(":")
                self._reverse[rid] = (rxn, detail.strip() or label.strip(), label.strip())

    def analyze(self, target: Molecule) -> list[RetroStep]:
        steps: list[RetroStep] = []
        seen: set[str] = set()
        for rid, (rxn, expl, name_fa) in self._reverse.items():
            try:
                product_sets = rxn.RunReactants((target.mol,))
            except Exception:
                continue
            for ps in product_sets:
                mols = []
                ok = True
                for p in ps:
                    try:
                        Chem.SanitizeMol(p)
                        mols.append(Molecule(p))
                    except Exception:
                        ok = False
                        break
                if not ok or not mols:
                    continue
                # skip trivial water-only fragments
                mols = [m for m in mols if m.formula != "H2O"]
                if not mols:
                    continue
                # dedupe on the precursor set itself (across all rules)
                key = ".".join(sorted(m.smiles for m in mols))
                if key in seen:
                    continue
                seen.add(key)
                steps.append(RetroStep(
                    rule_id=rid,
                    explanation_fa=f"{name_fa}: {expl}",
                    precursors=mols,
                ))
        return steps

    def base_building_blocks(self, target: Molecule, max_depth: int = 3) -> dict:
        """Recursively disconnect until no rule applies → 'base components'."""
        leaves: list[str] = []
        tree = self._recurse(target, max_depth, leaves, set())
        # unique leaves
        uniq = sorted(set(leaves))
        return {"tree": tree, "base_components": uniq}

    def _recurse(self, mol: Molecule, depth: int, leaves: list[str], visited: set[str]) -> dict:
        node = {"molecule": mol.formula, "smiles": mol.smiles, "children": []}
        if depth <= 0 or mol.smiles in visited:
            leaves.append(mol.smiles)
            return node
        visited = visited | {mol.smiles}
        steps = self.analyze(mol)
        if not steps:
            leaves.append(mol.smiles)
            return node
        # follow the first disconnection
        for p in steps[0].precursors:
            node["children"].append(self._recurse(p, depth - 1, leaves, visited))
        node["rule"] = steps[0].rule_id
        return node
