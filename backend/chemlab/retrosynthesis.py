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


def _clean_product(mol: Chem.Mol) -> Chem.Mol | None:
    """Sanitize a reaction product, repairing H bookkeeping if needed.

    Oxidation-state-changing transforms (alcohol⟶carbonyl, amine C–N cut) can
    leave a carbon with a stale explicit-H count and thus an over-valence. If
    the first sanitize fails, reset implicit-H handling and retry.
    """
    m = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(m)
        return m
    except Exception:
        pass
    m = Chem.Mol(mol)
    for atom in m.GetAtoms():
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    try:
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return None


# Reactions that form a new skeleton (worth reversing). Pure additions of small
# molecules (hydration etc.) are reversible too but produce trivial precursors.
_RETRO_RULES = {
    "esterification": "استر را به کربوکسیلیک‌اسید + الکل تجزیه می‌کند",
    "amide_formation": "آمید را به کربوکسیلیک‌اسید + آمین تجزیه می‌کند",
}

# Dedicated, more permissive retro disconnections (independent of forward rules).
# Using [#6] matches both aromatic and aliphatic carbons, so aryl esters
# (e.g. aspirin's acetate) are disconnected correctly.
#
# Two kinds of transforms live here:
#   * bond cleavages (2 precursors) — break the skeleton into building blocks
#   * functional-group interconversions / FGI (1 precursor) — change a group
#     without breaking the skeleton (e.g. nitro -> amine, alcohol -> carbonyl)
_GENERIC_DISCONNECTIONS = [
    # --- skeleton cleavages -------------------------------------------------
    ("ester_cut", "برش پیوند استری: استر ⟶ کربوکسیلیک‌اسید + الکل/فنل",
     "[C:1](=[O:2])[O:3][#6:4]>>[C:1](=[O:2])[OX2H:3].[#6:4][OX2H]"),
    ("amide_cut", "برش پیوند آمیدی: آمید ⟶ کربوکسیلیک‌اسید + آمین",
     "[C:1](=[O:2])[NX3:3]>>[C:1](=[O:2])[OX2H].[NX3:3][H]"),
    ("ether_cut", "برش پیوند اتری: اتر ⟶ دو الکل/فنل",
     "[#6:1][OX2;!$(O[CX3]=[OX1]):2][#6:3]>>[#6:1][OX2:2][H].[#6:3][OX2H]"),
    ("reductive_amination", "برش C–N آمین (احیای آمیناسیون): آمین ⟶ ترکیب کربونیل + آمین",
     "[CX4;H1,H2;!$([CX4][OX2]);!$([CX4]([#7])[#7]):1][NX3;!$([NX3][CX3]=[OX1]);!+:2]"
     ">>[CX3:1]=[O].[NX3:2][H]"),
    ("friedel_crafts_acyl", "برش C–C (آسیله‌شدن فریدل-کرافتس): آریل‌کتون ⟶ آرن + آسیل‌کلرید",
     "[c:1][CX3:2](=[OX1:3])[#6:4]>>[c:1][H].[Cl][C:2](=[O:3])[#6:4]"),
    ("aldol_cc", "برش C–C (آلدُل): بتا-هیدروکسی‌کربونیل ⟶ دو ترکیب کربونیل",
     "[CX4:1][CX4:2]([OX2H])[CX4:3][CX3:4]=[OX1:5]"
     ">>[CX3:1]=[O].[CX4:3][CX3:4]=[O:5]"),
    # --- functional-group interconversions (FGI) ---------------------------
    ("nitro_reduction", "احیای نیترو (معکوس): آمین آروماتیک ⟸ نیتروآرن",
     "[c:1][NX3;H1,H2;!$([NX3][CX3]=[OX1]):2]>>[c:1][N+:2](=O)[O-]"),
    ("alcohol_oxidation", "اکسایش (معکوس): الکل ⟸ ترکیب کربونیل",
     "[CX4;H1,H2;!$([CX4]([OX2,NX3])[OX2,NX3]);!$([CX4][OX2][CX3]=O):1][OX2H]"
     ">>[CX3:1]=[O]"),
]

# Transforms that merely break the molecule into smaller building blocks.
_CLEAVAGE_RULES = {"ester_cut", "amide_cut", "ether_cut", "reductive_amination",
                   "friedel_crafts_acyl", "aldol_cc"}


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
                    cleaned = _clean_product(p)
                    if cleaned is None:
                        ok = False
                        break
                    mols.append(Molecule(cleaned))
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
        # Only follow true skeleton cleavages (>1 fragment); FGIs don't shrink
        # the molecule and could otherwise loop. They still appear in analyze().
        cleavages = [s for s in steps
                     if s.rule_id in _CLEAVAGE_RULES and len(s.precursors) > 1]
        if not cleavages:
            leaves.append(mol.smiles)
            return node
        chosen = cleavages[0]
        for p in chosen.precursors:
            node["children"].append(self._recurse(p, depth - 1, leaves, visited))
        node["rule"] = chosen.rule_id
        return node
