"""
nlp.py — Natural-language scenario parsing (Persian + English).

Turns a free-text scenario like
    «استیک اسید و اتانول را با کاتالیزور سولفوریک اسید در دمای ۸۰ درجه حرارت بده»
into structured reactants + Conditions for the reaction engine. Heuristic,
dictionary-driven: it scans the text for known compound names, numeric
temperatures, catalysts, light and technique keywords.
"""
from __future__ import annotations

import re

from .conditions import Conditions, Technique
from .known_compounds import _INDEX

# Persian digit normalisation.
_PERSIAN_DIGITS = str.maketrans("۰۱۲۳۴۵۶۷۸۹٫", "0123456789.")

_TECH_KEYWORDS = {
    "تقطیر": Technique.DISTILLATION, "distill": Technique.DISTILLATION,
    "رفلاکس": Technique.REFLUX, "بازروان": Technique.REFLUX, "reflux": Technique.REFLUX,
    "برق‌کافت": Technique.ELECTROLYSIS, "الکترولیز": Technique.ELECTROLYSIS,
    "electrolys": Technique.ELECTROLYSIS,
    "نورکافت": Technique.PHOTOLYSIS, "photolys": Technique.PHOTOLYSIS,
    "تبلور": Technique.CRYSTALLIZATION, "crystall": Technique.CRYSTALLIZATION,
    "صاف": Technique.FILTRATION, "filtrat": Technique.FILTRATION,
    "حرارت": Technique.HEATING, "گرم": Technique.HEATING, "heat": Technique.HEATING,
    "سرد": Technique.COOLING, "cool": Technique.COOLING,
}

_CATALYST_HINTS = {
    "سولفوریک": "H2SO4", "h2so4": "H2SO4", "sulfuric": "H2SO4",
    "کلریدریک": "HCl", "hcl": "HCl",
    "فسفریک": "H3PO4", "h3po4": "H3PO4",
    "پالادیم": "Pd", "pd": "Pd", "پلاتین": "Pt", "نیکل": "Ni", "ni": "Ni",
    "سود": "NaOH", "naoh": "NaOH", "هیدروکسید سدیم": "NaOH",
    "alcl3": "AlCl3", "کلرید آلومینیوم": "AlCl3",
}


def _normalise(text: str) -> str:
    return text.translate(_PERSIAN_DIGITS)


def parse_scenario(text: str) -> dict:
    raw = text or ""
    norm = _normalise(raw).lower()

    # --- reactants: longest known names first to avoid partial matches ---
    reactants: list[str] = []
    found_spans: list[tuple[int, int]] = []
    for name in sorted(_INDEX.keys(), key=len, reverse=True):
        if len(name) < 2:
            continue
        idx = norm.find(name)
        if idx == -1:
            continue
        # skip if overlapping an already-matched span
        if any(s <= idx < e or s < idx + len(name) <= e for s, e in found_spans):
            continue
        # skip if it directly follows a "catalyst" cue → it's the catalyst, not a reactant
        prefix = norm[max(0, idx - 14):idx]
        if any(cue in prefix for cue in ("کاتالیزور", "catalyst", "کاتالیست")):
            found_spans.append((idx, idx + len(name)))
            continue
        found_spans.append((idx, idx + len(name)))
        reactants.append(name)

    # --- temperature: number followed by درجه/°C/c ---
    temp = 25.0
    mt = re.search(r"(-?\d+(?:\.\d+)?)\s*(?:درجه|°c|°|c\b|سلسیوس|celsius)", norm)
    if mt:
        temp = float(mt.group(1))

    # --- catalyst ---
    catalyst = None
    for hint, cat in _CATALYST_HINTS.items():
        if hint in norm:
            catalyst = cat
            break

    # --- light & techniques ---
    light = any(w in norm for w in ("نور", "uv", "فرابنفش", "light"))
    techs: list[Technique] = []
    for kw, tech in _TECH_KEYWORDS.items():
        if kw in norm and tech not in techs:
            techs.append(tech)
    if any(t in techs for t in (Technique.HEATING,)) and temp == 25.0:
        temp = 80.0  # "حرارت بده" implies heating

    cond = Conditions(temperature_c=temp, catalyst=catalyst,
                      light=light, techniques=techs)
    return {
        "reactants": reactants,
        "conditions": cond.to_dict(),
        "_conditions_obj": cond,
        "parsed_fa": (
            f"مواد شناسایی‌شده: {('، '.join(reactants)) or 'هیچ'}؛ "
            f"دما: {temp}°C؛ کاتالیزور: {catalyst or 'ندارد'}؛ "
            f"نور: {'بله' if light else 'خیر'}؛ "
            f"تکنیک‌ها: {('، '.join(t.value for t in techs)) or 'ندارد'}."
        ),
    }
