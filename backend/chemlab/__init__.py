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
from .known_compounds import all_compounds, resolve_name, all_drugs
from .balance import balance
from .thermo import thermodynamics, kinetics
from .stoichiometry import stoichiometry
from .solution import ph_of_solution, solubility, lookup_pka
from .separation import distill, boiling_point
from .redox import oxidation_states, cell_potential, electrolysis, STANDARD_POTENTIALS
from .admet import admet, interaction
from .spectra import full_spectra, ir_spectrum, nmr_h_spectrum, mass_spectrum
from .nlp import parse_scenario

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
    "all_compounds", "resolve_name", "all_drugs",
    "balance", "thermodynamics", "kinetics", "stoichiometry",
    "ph_of_solution", "solubility", "lookup_pka",
    "distill", "boiling_point",
    "oxidation_states", "cell_potential", "electrolysis", "STANDARD_POTENTIALS",
    "admet", "interaction",
    "full_spectra", "ir_spectrum", "nmr_h_spectrum", "mass_spectrum",
    "parse_scenario",
]
