"""
conditions.py — Reaction / environment conditions.

Captures the physical environment of an experiment: temperature, pressure,
humidity, presence of catalyst, solvent, light, time, and special techniques
(distillation, reflux, electrolysis, ...). These conditions gate which
reactions in the engine are feasible and shift outcomes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Technique(str, Enum):
    NONE = "none"
    HEATING = "heating"          # گرمادهی
    COOLING = "cooling"          # سرمادهی
    DISTILLATION = "distillation"  # تقطیر
    REFLUX = "reflux"            # بازروانی
    ELECTROLYSIS = "electrolysis"  # برق‌کافت
    PHOTOLYSIS = "photolysis"    # نورکافت
    CRYSTALLIZATION = "crystallization"  # تبلور
    FILTRATION = "filtration"    # صاف‌کردن


TECHNIQUE_FA = {
    Technique.NONE: "بدون تکنیک خاص",
    Technique.HEATING: "گرمادهی",
    Technique.COOLING: "سرمادهی",
    Technique.DISTILLATION: "تقطیر",
    Technique.REFLUX: "بازروانی (رفلاکس)",
    Technique.ELECTROLYSIS: "برق‌کافت",
    Technique.PHOTOLYSIS: "نورکافت",
    Technique.CRYSTALLIZATION: "تبلور",
    Technique.FILTRATION: "صاف‌کردن",
}


@dataclass
class Conditions:
    """Environmental conditions of an experiment."""

    temperature_c: float = 25.0          # دما (سلسیوس)
    pressure_atm: float = 1.0            # فشار (اتمسفر)
    humidity_pct: float = 50.0           # رطوبت نسبی ٪
    ph: Optional[float] = None           # pH محیط
    catalyst: Optional[str] = None       # کاتالیزور
    solvent: Optional[str] = None        # حلال
    light: bool = False                  # نور (UV/مرئی)
    duration_min: Optional[float] = None  # مدت زمان (دقیقه)
    techniques: list[Technique] = field(default_factory=list)

    @property
    def temperature_k(self) -> float:
        return self.temperature_c + 273.15

    def has(self, tech: Technique) -> bool:
        return tech in self.techniques

    @property
    def is_heated(self) -> bool:
        return self.temperature_c >= 60 or self.has(Technique.HEATING) or self.has(Technique.REFLUX)

    @property
    def is_cold(self) -> bool:
        return self.temperature_c <= 5 or self.has(Technique.COOLING)

    def describe_fa(self) -> str:
        parts = [f"دما {self.temperature_c}°C", f"فشار {self.pressure_atm} atm"]
        if self.ph is not None:
            parts.append(f"pH={self.ph}")
        if self.catalyst:
            parts.append(f"کاتالیزور: {self.catalyst}")
        if self.solvent:
            parts.append(f"حلال: {self.solvent}")
        if self.light:
            parts.append("در حضور نور")
        for t in self.techniques:
            parts.append(TECHNIQUE_FA.get(t, t.value))
        return "، ".join(parts)

    def to_dict(self) -> dict:
        return {
            "temperature_c": self.temperature_c,
            "temperature_k": round(self.temperature_k, 2),
            "pressure_atm": self.pressure_atm,
            "humidity_pct": self.humidity_pct,
            "ph": self.ph,
            "catalyst": self.catalyst,
            "solvent": self.solvent,
            "light": self.light,
            "duration_min": self.duration_min,
            "techniques": [t.value for t in self.techniques],
            "is_heated": self.is_heated,
            "is_cold": self.is_cold,
            "description_fa": self.describe_fa(),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Conditions":
        techs = []
        for t in d.get("techniques", []) or []:
            try:
                techs.append(Technique(t))
            except ValueError:
                pass
        return cls(
            temperature_c=float(d.get("temperature_c", 25.0)),
            pressure_atm=float(d.get("pressure_atm", 1.0)),
            humidity_pct=float(d.get("humidity_pct", 50.0)),
            ph=(float(d["ph"]) if d.get("ph") not in (None, "") else None),
            catalyst=d.get("catalyst") or None,
            solvent=d.get("solvent") or None,
            light=bool(d.get("light", False)),
            duration_min=(float(d["duration_min"]) if d.get("duration_min") not in (None, "") else None),
            techniques=techs,
        )
