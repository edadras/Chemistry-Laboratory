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


def _needs_lewis_acid(c: Conditions) -> tuple[bool, str]:
    cat = (c.catalyst or "").lower()
    if any(k in cat for k in ["fecl3", "alcl3", "febr3", "لوویس", "آهن", "آلومینیوم"]):
        return True, "کاتالیزور اسید لوویس فراهم است"
    return False, "این واکنش به کاتالیزور اسید لوویس (FeCl3/AlCl3) نیاز دارد"


def _needs_acid_and_heat(c: Conditions) -> tuple[bool, str]:
    return _needs_acid_catalyst(c)


# NOTE: SMARTS are validated at import-time by reactions.py.
REACTION_RULES: list[ReactionRule] = [
    ReactionRule(
        rid="esterification",
        name="Fischer esterification / O-acylation",
        name_fa="استری‌شدن (الکل و فنل)",
        smarts="[C:1](=[O:2])[OX2H].[OX2H:4][#6:5]>>[C:1](=[O:2])[O:4][#6:5].[OH2]",
        category="organic",
        description_fa="کربوکسیلیک اسید + الکل/فنل ⟶ استر + آب (با کاتالیزور اسیدی و گرما)",
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
    ReactionRule(
        rid="hydrohalogenation",
        name="Hydrohalogenation (Markovnikov)",
        name_fa="افزایش هیدروهالید به آلکن",
        smarts="[C:1]=[C:2].[Cl;H1:3]>>[C:1][C:2][Cl:3]",
        category="organic",
        description_fa="آلکن + HCl ⟶ آلکیل‌هالید (افزایش مارکوونیکوف)",
        feasible=_always,
    ),
    ReactionRule(
        rid="dehydration",
        name="Acid-catalysed dehydration",
        name_fa="آب‌زدایی الکل (تشکیل آلکن)",
        smarts="[CX4;H1,H2:1][CX4:2][OX2H]>>[CX3:1]=[CX3:2].[OH2]",
        category="organic",
        description_fa="الکل ⟶ آلکن + آب (کاتالیزور اسیدی و گرما)",
        feasible=_needs_acid_catalyst,
        needs_catalyst="H2SO4",
    ),
    ReactionRule(
        rid="carbonyl_reduction",
        name="Carbonyl reduction",
        name_fa="احیای کربونیل به الکل",
        smarts="[CX3:1]=[OX1:2].[H][H]>>[CX4:1][OX2H:2]",
        category="organic",
        description_fa="آلدهید/کتون + H2 ⟶ الکل (کاتالیزور فلزی یا NaBH4)",
        feasible=_needs_metal_catalyst,
        needs_catalyst="Ni/NaBH4",
    ),
    ReactionRule(
        rid="aromatic_nitration",
        name="Aromatic nitration",
        name_fa="نیترودارشدن حلقه‌ی آروماتیک",
        smarts="[cH:1].O[N+](=O)[O-]>>[c:1][N+](=O)[O-].[OH2]",
        category="organic",
        description_fa="بنزن + HNO3 (با H2SO4) ⟶ نیتروبنزن + آب",
        feasible=_needs_acid_catalyst,
        needs_catalyst="H2SO4",
    ),
    ReactionRule(
        rid="aromatic_halogenation",
        name="Electrophilic aromatic halogenation",
        name_fa="هالوژن‌دارشدن آروماتیک (الکتروفیلی)",
        smarts="[cH:1].[Cl][Cl]>>[c:1][Cl].[Cl][H]",
        category="organic",
        description_fa="بنزن + Cl2 (با FeCl3) ⟶ کلروبنزن + HCl",
        feasible=_needs_lewis_acid,
        needs_catalyst="FeCl3",
    ),
    ReactionRule(
        rid="ester_hydrolysis_acid",
        name="Acidic ester hydrolysis",
        name_fa="هیدرولیز اسیدی استر",
        smarts="[C:1](=[O:2])[O:3][#6:4].[OH2]>>[C:1](=[O:2])[O:3][H].[#6:4][OX2H]",
        category="organic",
        description_fa="استر + آب (محیط اسیدی) ⟶ کربوکسیلیک‌اسید + الکل",
        feasible=_needs_acid_catalyst,
    ),
    ReactionRule(
        rid="decarboxylation",
        name="Decarboxylation",
        name_fa="کربوکسیل‌زدایی",
        smarts="[#6:1][CX3](=O)[OX2H]>>[#6:1][H].O=C=O",
        category="organic",
        description_fa="کربوکسیلیک‌اسید ⟶ هیدروکربن + CO2 (با گرما)",
        feasible=_needs_heat,
    ),
]
