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


def test_retro_reductive_amination():
    retro = Retrosynthesizer()
    phenylephrine = Molecule.from_smiles("CNC[C@H](O)c1cccc(O)c1")
    steps = retro.analyze(phenylephrine)
    rule_ids = {s.rule_id for s in steps}
    assert "reductive_amination" in rule_ids
    blocks = retro.base_building_blocks(phenylephrine)["base_components"]
    assert len(blocks) >= 2  # broke into building blocks (not just itself)


def test_synthesize_paracetamol_confirmed():
    sp = ScenarioProcessor()
    res = sp.synthesize("استامینوفن")
    assert res["ok"]
    assert any(r["forward_confirmed"] for r in res["routes"])
    assert res["drug_info"]["category"] == "analgesic"


def test_synthesize_aspirin_confirmed():
    sp = ScenarioProcessor()
    res = sp.synthesize("آسپرین")
    assert any(r["disconnection"] == "ester_cut" and r["forward_confirmed"]
               for r in res["routes"])


def test_drug_database_metadata():
    from chemlab.known_compounds import all_drugs, drug_info
    drugs = all_drugs()
    assert len(drugs) >= 25
    cats = {d["category"] for d in drugs}
    assert {"antibiotic", "antiviral", "analgesic"} <= cats
    info = drug_info("ciprofloxacin")
    assert info and info["category"] == "antibiotic"


def test_balance_combustion():
    from chemlab.balance import balance
    res = balance([Molecule.parse("propane"), Molecule.from_smiles("O=O")],
                  [Molecule.from_smiles("O=C=O"), Molecule.from_smiles("O")])
    assert res is not None
    assert res["reactant_coeffs"] == [1, 5]
    assert res["product_coeffs"] == [3, 4]


def test_balance_water():
    from chemlab.balance import balance
    res = balance([Molecule.from_smiles("[H][H]"), Molecule.from_smiles("O=O")],
                  [Molecule.from_smiles("O")])
    assert res["reactant_coeffs"] == [2, 1] and res["product_coeffs"] == [2]


def test_outcome_has_balanced_equation():
    eng = ReactionEngine()
    out = eng.predict([Molecule.parse("methane"), Molecule.from_smiles("O=O")],
                      Conditions(temperature_c=600))
    comb = [o for o in out if o.rule_id == "combustion"][0]
    d = comb.to_dict(include_svg=False)
    assert d["balanced_equation"] == "CH4 + 2 O2 → CO2 + 2 H2O"


def test_multistep_synthesis():
    sp = ScenarioProcessor()
    # phenacetin decomposes in two steps (ether synthesis + amidation)
    res = sp.synthesize("CCOc1ccc(NC(C)=O)cc1", multistep=True, max_steps=4)
    assert res["ok"] and res["multistep"]
    assert len(res["steps"]) >= 2
    assert any(s["forward_confirmed"] for s in res["steps"])
    assert res["starting_materials"]


def test_expanded_drug_count():
    from chemlab.known_compounds import all_drugs
    assert len(all_drugs()) >= 40


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
