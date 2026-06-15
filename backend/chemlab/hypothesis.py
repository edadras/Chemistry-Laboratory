"""
hypothesis.py — Autonomous hypothesis generation & testing.

The engine forms hypotheses ("if I combine A and B under these conditions, X
should form") and then *tests* each one by running the forward reaction engine,
evaluating the product, and scoring it. In pharma mode products are scored by
drug-likeness (QED) so the system can propose promising candidates on its own.

This is a rule-based, transparent surrogate for a real experiment loop: every
hypothesis carries its prediction, the test result, and a verdict.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

from .conditions import Conditions
from .molecule import Molecule
from .reactions import ReactionEngine, Outcome
from .pharma import DrugProfile


@dataclass
class Hypothesis:
    statement_fa: str
    reactants: list[str]          # formulas
    predicted_products: list[str]  # formulas
    feasible: bool
    score: float                  # 0..1 desirability (QED-based in pharma mode)
    verdict_fa: str
    outcome: Outcome = field(repr=False, default=None)

    def to_dict(self, include_svg: bool = False) -> dict:
        d = {
            "statement_fa": self.statement_fa,
            "reactants": self.reactants,
            "predicted_products": self.predicted_products,
            "feasible": self.feasible,
            "score": round(self.score, 3),
            "verdict_fa": self.verdict_fa,
        }
        if self.outcome is not None:
            d["outcome"] = self.outcome.to_dict(include_svg=include_svg)
        return d


class HypothesisEngine:
    def __init__(self) -> None:
        self.engine = ReactionEngine()

    def generate(self, reactants: list[Molecule], cond: Conditions,
                 pharma_mode: bool = False, top_k: int = 10) -> list[Hypothesis]:
        """Enumerate pairwise/triple combinations and test each hypothesis."""
        hyps: list[Hypothesis] = []
        seen: set[str] = set()
        n = len(reactants)
        sizes = [s for s in (2, 3) if s <= n] or [n]
        for size in sizes:
            for combo in itertools.combinations(reactants, size):
                outcomes = self.engine.predict(list(combo), cond, include_infeasible=True)
                for oc in outcomes:
                    key = oc.rule_id + "|" + "+".join(sorted(m.formula for m in combo))
                    if key in seen:
                        continue
                    seen.add(key)
                    hyps.append(self._evaluate(combo, oc, pharma_mode))
        # rank: feasible first, then by score
        hyps.sort(key=lambda h: (not h.feasible, -h.score))
        return hyps[:top_k]

    def _evaluate(self, combo, oc: Outcome, pharma_mode: bool) -> Hypothesis:
        r_formulas = [m.formula for m in combo]
        p_formulas = [p.formula for p in oc.products]
        score = 0.5
        verdict = oc.reason_fa
        if pharma_mode and oc.products:
            # score by best product drug-likeness
            best = max(oc.products, key=lambda p: _safe_qed(p))
            q = _safe_qed(best)
            score = q
            prof = DrugProfile(best)
            verdict = (f"بهترین محصول {best.formula} با QED={q:.2f}؛ "
                       f"{prof.verdict_fa()}")
        else:
            score = 0.7 if oc.feasible else 0.3
        statement = (f"فرضیه: ترکیب {' + '.join(r_formulas)} تحت شرایط داده‌شده "
                     f"از طریق «{oc.rule_name_fa}» به {' + '.join(p_formulas)} منجر می‌شود.")
        return Hypothesis(
            statement_fa=statement,
            reactants=r_formulas,
            predicted_products=p_formulas,
            feasible=oc.feasible,
            score=score,
            verdict_fa=verdict,
            outcome=oc,
        )


def _safe_qed(mol: Molecule) -> float:
    try:
        from rdkit.Chem import QED
        return float(QED.qed(mol.mol))
    except Exception:
        return 0.0
