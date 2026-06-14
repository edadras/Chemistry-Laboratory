"""Smoke + behavioural tests for the chemlab engine."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from chemlab import (
    periodic_table, Molecule, Conditions, ReactionEngine,
    Retrosynthesizer, ScenarioProcessor, DrugProfile,
)
from chemlab.conditions import Technique


def test_periodic_table():
    fe = periodic_table.by_symbol("Fe")
    assert fe.atomic_number == 26
    assert fe.name == "Iron"
    assert periodic_table.by_number(8).symbol == "O"
    assert len(periodic_table.all()) == 118


def test_molecule_parse_and_descriptors():
    asp = Molecule.parse("آسپرین")
    assert asp.formula == "C9H8O4"
    d = asp.descriptors()
    assert d["molar_mass"] > 150
    assert "svg" in asp.to_dict()


def test_esterification():
    eng = ReactionEngine()
    acid = Molecule.parse("acetic acid")
    etoh = Molecule.parse("ethanol")
    cond = Conditions(temperature_c=80, catalyst="H2SO4")
    outcomes = eng.predict([acid, etoh], cond)
    ester = [o for o in outcomes if o.rule_id == "esterification"]
    assert ester and ester[0].feasible
    formulas = [p.formula for p in ester[0].products]
    assert "C4H8O2" in formulas  # ethyl acetate


def test_esterification_needs_conditions():
    eng = ReactionEngine()
    acid = Molecule.parse("acetic acid")
    etoh = Molecule.parse("ethanol")
    outcomes = eng.predict([acid, etoh], Conditions())  # no catalyst/heat
    ester = [o for o in outcomes if o.rule_id == "esterification"]
    assert ester and not ester[0].feasible


def test_neutralization():
    eng = ReactionEngine()
    hcl = Molecule.parse("hydrochloric acid")
    naoh = Molecule.parse("sodium hydroxide")
    outcomes = eng.predict([hcl, naoh], Conditions())
    neut = [o for o in outcomes if o.rule_id == "neutralization"]
    assert neut and neut[0].feasible


def test_combustion():
    eng = ReactionEngine()
    ch4 = Molecule.parse("methane")
    o2 = Molecule.from_smiles("O=O")
    outcomes = eng.predict([ch4, o2], Conditions(temperature_c=600))
    comb = [o for o in outcomes if o.rule_id == "combustion"]
    assert comb
    formulas = {p.formula for p in comb[0].products}
    assert {"CO2", "H2O"} <= formulas


def test_retrosynthesis_aspirin():
    retro = Retrosynthesizer()
    asp = Molecule.parse("آسپرین")
    steps = retro.analyze(asp)
    assert steps
    # acetic acid + salicylic acid should appear among precursors
    all_formulas = {p.formula for s in steps for p in s.precursors}
    assert "C2H4O2" in all_formulas and "C7H6O3" in all_formulas


def test_hypothesis_discovers_aspirin():
    sp = ScenarioProcessor()
    res = sp.hypothesize(
        ["acetic acid", "salicylic acid"],
        Conditions(temperature_c=80, catalyst="H2SO4"),
        pharma_mode=True, top_k=5,
    )
    assert res["ok"]
    products = {f for h in res["hypotheses"] for f in h["predicted_products"]}
    assert "C9H8O4" in products  # aspirin discovered


def test_drug_profile():
    asp = Molecule.parse("آسپرین")
    prof = DrugProfile(asp).to_dict()
    assert prof["lipinski"]["passes"] is True
    assert 0 <= prof["qed"] <= 1


def test_scenario_react_summary():
    sp = ScenarioProcessor()
    res = sp.react(["acetic acid", "ethanol"],
                   Conditions(temperature_c=80, catalyst="H2SO4"))
    assert res["ok"] and res["main_product"]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {fn.__name__}: {e}")
        except Exception as e:
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(fns)} tests passed")
    sys.exit(0 if passed == len(fns) else 1)
