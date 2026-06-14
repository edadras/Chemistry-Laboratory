"""
chemlab — a computational chemistry laboratory engine.

Public surface:
    PeriodicTable / periodic_table  — elements & their properties
    Molecule                        — structure, bonds, descriptors, depiction
    Conditions / Technique          — experiment environment
    ReactionEngine                  — forward reaction prediction
    Retrosynthesizer                — reverse engineering of a target
    HypothesisEngine                — autonomous hypothesis generation/testing
    DrugProfile                     — pharmaceutical analysis
    ScenarioProcessor               — high-level orchestration for the API
"""
from .elements import PeriodicTable, periodic_table, Element
from .molecule import Molecule, MoleculeError
from .conditions import Conditions, Technique
from .reactions import ReactionEngine, Outcome
from .retrosynthesis import Retrosynthesizer
from .hypothesis import HypothesisEngine, Hypothesis
from .pharma import DrugProfile, indication_for
from .scenario import ScenarioProcessor
from .known_compounds import all_compounds, resolve_name

__version__ = "0.1.0"

__all__ = [
    "PeriodicTable", "periodic_table", "Element",
    "Molecule", "MoleculeError",
    "Conditions", "Technique",
    "ReactionEngine", "Outcome",
    "Retrosynthesizer",
    "HypothesisEngine", "Hypothesis",
    "DrugProfile", "indication_for",
    "ScenarioProcessor",
    "all_compounds", "resolve_name",
]
