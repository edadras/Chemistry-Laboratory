"""
reaction_data.py — Curated library of named reaction rules.

Each rule is an RDKit reaction-SMARTS transform plus the conditions under
which it is feasible. The `feasible` callable receives a `Conditions` object
and returns (ok, reason_fa). The engine uses these to do forward prediction
and (by analysing precursors) retrosynthesis.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

from .conditions import Conditions, Technique


@dataclass
class ReactionRule:
    rid: str
    name: str
    name_fa: str
    smarts: str                 # RDKit reaction SMARTS (a>>b)
    category: str               # organic / inorganic / ...
    description_fa: str
    feasible: Callable[[Conditions], tuple[bool, str]]
    needs_catalyst: Optional[str] = None


def _always(_: Conditions) -> tuple[bool, str]:
    return True, "در شرایط استاندارد انجام می‌شود"


def _needs_heat(c: Conditions) -> tuple[bool, str]:
    if c.is_heated:
        return True, "گرما فراهم است"
    return False, "این واکنش به گرمادهی نیاز دارد (دمای بالاتر یا تکنیک گرمادهی)"


def _needs_acid_catalyst(c: Conditions) -> tuple[bool, str]:
    cat = (c.catalyst or "").lower()
    acidic = any(k in cat for k in ["h2so4", "acid", "اسید", "hcl", "h3po4"])
    if acidic or (c.ph is not None and c.ph < 3):
        if c.is_heated:
            return True, "کاتالیزور اسیدی و گرما فراهم است"
        return False, "کاتالیزور اسیدی هست ولی واکنش به گرما هم نیاز دارد"
    return False, "این واکنش به کاتالیزور اسیدی (مثل H2SO4) و گرما نیاز دارد"


def _needs_metal_catalyst(c: Conditions) -> tuple[bool, str]:
    cat = (c.catalyst or "").lower()
    if any(k in cat for k in ["pd", "pt", "ni", "پالادیم", "پلاتین", "نیکل"]):
        return True, "کاتالیزور فلزی فراهم است"
    return False, "هیدروژناسیون به کاتالیزور فلزی (Pd/Pt/Ni) نیاز دارد"


def _needs_base(c: Conditions) -> tuple[bool, str]:
    if c.ph is not None and c.ph > 9:
        return True, "محیط بازی فراهم است"
    cat = (c.catalyst or "").lower()
    if any(k in cat for k in ["naoh", "koh", "oh-", "باز", "سود"]):
        return True, "باز قوی فراهم است"
    return False, "این واکنش به محیط بازی (مثل NaOH) نیاز دارد"


def _needs_light(c: Conditions) -> tuple[bool, str]:
    if c.light or c.has(Technique.PHOTOLYSIS):
        return True, "نور فراهم است"
    return False, "این واکنش رادیکالی به نور (UV) نیاز دارد"


# NOTE: SMARTS are validated at import-time by reactions.py.
REACTION_RULES: list[ReactionRule] = [
    ReactionRule(
        rid="esterification",
        name="Fischer esterification",
        name_fa="استری‌شدن فیشر",
        smarts="[C:1](=[O:2])[OX2H].[OX2H:4][C:5]>>[C:1](=[O:2])[O:4][C:5].[OH2]",
        category="organic",
        description_fa="کربوکسیلیک اسید + الکل ⟶ استر + آب (با کاتالیزور اسیدی و گرما)",
        feasible=_needs_acid_catalyst,
        needs_catalyst="H2SO4",
    ),
    ReactionRule(
        rid="amide_formation",
        name="Amide formation",
        name_fa="تشکیل آمید",
        smarts="[C:1](=[O:2])[OX2H].[NX3;H2,H1:4]>>[C:1](=[O:2])[N:4].[OH2]",
        category="organic",
        description_fa="کربوکسیلیک اسید + آمین ⟶ آمید + آب (با گرما)",
        feasible=_needs_heat,
    ),
    ReactionRule(
        rid="hydrogenation",
        name="Catalytic hydrogenation",
        name_fa="هیدروژناسیون کاتالیزوری",
        smarts="[C:1]=[C:2].[H][H]>>[C:1][C:2]",
        category="organic",
        description_fa="آلکن + H2 ⟶ آلکان (با کاتالیزور فلزی Pd/Pt/Ni)",
        feasible=_needs_metal_catalyst,
        needs_catalyst="Pd/Pt/Ni",
    ),
    ReactionRule(
        rid="alkene_hydration",
        name="Acid-catalyzed hydration",
        name_fa="آب‌گیری آلکن (هیدراسیون)",
        smarts="[C:1]=[C:2].[OH2]>>[C:1][C:2][OX2H]",
        category="organic",
        description_fa="آلکن + آب ⟶ الکل (با کاتالیزور اسیدی)",
        feasible=_needs_acid_catalyst,
        needs_catalyst="H3PO4/H2SO4",
    ),
    ReactionRule(
        rid="alkene_halogenation",
        name="Halogen addition",
        name_fa="افزایش هالوژن به آلکن",
        smarts="[C:1]=[C:2].[Cl:3][Cl:4]>>[C:1]([Cl:3])[C:2][Cl:4]",
        category="organic",
        description_fa="آلکن + Cl2 ⟶ دی‌هالید (افزایشی، سریع)",
        feasible=_always,
    ),
    ReactionRule(
        rid="saponification",
        name="Saponification",
        name_fa="صابونی‌شدن (هیدرولیز بازی استر)",
        smarts="[C:1](=[O:2])[O:3][C:4].[OH2]>>[C:1](=[O:2])[O:3].[OX2H][C:4]",
        category="organic",
        description_fa="استر + آب (محیط بازی) ⟶ کربوکسیلات + الکل",
        feasible=_needs_base,
    ),
    ReactionRule(
        rid="alcohol_oxidation",
        name="Oxidation of primary alcohol",
        name_fa="اکسایش الکل نوع اول",
        smarts="[CX4;H2:1][OX2H]>>[CX3:1]=[O]",
        category="organic",
        description_fa="الکل نوع اول ⟶ آلدهید (اکسنده + گرما)",
        feasible=_needs_heat,
    ),
    ReactionRule(
        rid="radical_halogenation",
        name="Free-radical halogenation",
        name_fa="هالوژن‌دارشدن رادیکالی",
        smarts="[CX4;H1,H2,H3:1].[Cl][Cl]>>[C:1][Cl]",
        category="organic",
        description_fa="آلکان + Cl2 (در حضور نور) ⟶ آلکیل‌هالید + HCl",
        feasible=_needs_light,
    ),
]
