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
from rdkit import RDLogger as _RDLogger

# RDKit's C++ layer is chatty (reaction atom-mapping notes, transient
# sanitization errors we already handle via try/except). Silence it so server
# and CLI output stays readable; failures still surface through Python.
_RDLogger.DisableLog("rdApp.*")

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
from .qsar import QSARModel, demo_model, featurize
from .generator import evolve, mutate, crossover
from .discovery import discover, candidate_profile, synth_accessibility
from .docking import dock
from .targets import list_targets, resolve_target
from . import workbench
from . import external_db

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
    "QSARModel", "demo_model", "featurize",
    "evolve", "mutate", "crossover",
    "discover", "candidate_profile", "synth_accessibility",
    "dock", "list_targets", "resolve_target", "workbench", "external_db",
]
