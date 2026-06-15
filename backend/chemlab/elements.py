"""
elements.py — Periodic table layer.

Provides rich access to the chemical elements. Base atomic data (symbol,
atomic number, mass, default valence, electronegativity, etc.) comes from
RDKit's built-in periodic table, augmented with a curated enrichment table
for physical/chemical properties that RDKit does not carry (melting/boiling
points, electron configuration, category, common oxidation states, discovery,
biological role, ...).

The goal is to be a *single source of truth* about elements for the rest of
the laboratory engine.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Optional

from rdkit import Chem

from .enrichment_data import ELEMENT_ENRICHMENT, CATEGORY_FA


@dataclass(frozen=True)
class Element:
    """A single chemical element with merged RDKit + curated data."""

    atomic_number: int
    symbol: str
    name: str
    atomic_mass: float
    # Curated / enriched fields (may be None when unknown).
    category: Optional[str] = None
    group: Optional[int] = None
    period: Optional[int] = None
    electron_config: Optional[str] = None
    electronegativity: Optional[float] = None
    melting_point_k: Optional[float] = None
    boiling_point_k: Optional[float] = None
    density: Optional[float] = None  # g/cm^3 (solids/liquids) or g/L (gases)
    oxidation_states: tuple[int, ...] = field(default_factory=tuple)
    common_valence: tuple[int, ...] = field(default_factory=tuple)
    phase_stp: Optional[str] = None  # solid / liquid / gas
    discovered: Optional[str] = None
    biological_role: Optional[str] = None

    @property
    def category_fa(self) -> str:
        """Persian label for the element category."""
        return CATEGORY_FA.get(self.category or "", self.category or "نامشخص")

    def to_dict(self) -> dict:
        return {
            "atomic_number": self.atomic_number,
            "symbol": self.symbol,
            "name": self.name,
            "atomic_mass": round(self.atomic_mass, 4),
            "category": self.category,
            "category_fa": self.category_fa,
            "group": self.group,
            "period": self.period,
            "electron_config": self.electron_config,
            "electronegativity": self.electronegativity,
            "melting_point_k": self.melting_point_k,
            "boiling_point_k": self.boiling_point_k,
            "density": self.density,
            "oxidation_states": list(self.oxidation_states),
            "common_valence": list(self.common_valence),
            "phase_stp": self.phase_stp,
            "discovered": self.discovered,
            "biological_role": self.biological_role,
        }


class PeriodicTable:
    """Façade over RDKit's periodic table + curated enrichment."""

    def __init__(self) -> None:
        self._rd = Chem.GetPeriodicTable()

    @lru_cache(maxsize=200)
    def by_number(self, z: int) -> Element:
        if not (1 <= z <= 118):
            raise ValueError(f"Atomic number out of range: {z}")
        symbol = self._rd.GetElementSymbol(z)
        return self._build(z, symbol)

    def by_symbol(self, symbol: str) -> Element:
        symbol = symbol.strip().capitalize()
        z = self._rd.GetAtomicNumber(symbol)
        if z == 0:
            raise ValueError(f"Unknown element symbol: {symbol!r}")
        return self.by_number(z)

    def by_name(self, name: str) -> Element:
        key = name.strip().lower()
        for z in range(1, 119):
            enr = ELEMENT_ENRICHMENT.get(z, {})
            if enr.get("name", "").lower() == key or enr.get("name_fa") == name.strip():
                return self.by_number(z)
        raise ValueError(f"Unknown element name: {name!r}")

    def _build(self, z: int, symbol: str) -> Element:
        enr = ELEMENT_ENRICHMENT.get(z, {})
        # Default valence list straight from RDKit.
        try:
            valences = tuple(v for v in self._rd.GetValenceList(z) if v > 0)
        except Exception:
            valences = ()
        return Element(
            atomic_number=z,
            symbol=symbol,
            name=enr.get("name", symbol),
            atomic_mass=self._rd.GetAtomicWeight(z),
            category=enr.get("category"),
            group=enr.get("group"),
            period=enr.get("period"),
            electron_config=enr.get("electron_config"),
            electronegativity=enr.get("electronegativity"),
            melting_point_k=enr.get("melting_point_k"),
            boiling_point_k=enr.get("boiling_point_k"),
            density=enr.get("density"),
            oxidation_states=tuple(enr.get("oxidation_states", ())),
            common_valence=valences,
            phase_stp=enr.get("phase_stp"),
            discovered=enr.get("discovered"),
            biological_role=enr.get("biological_role"),
        )

    def all(self) -> list[Element]:
        return [self.by_number(z) for z in range(1, 119)]


# Singleton convenience instance.
periodic_table = PeriodicTable()
