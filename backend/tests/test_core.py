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


def test_thermo_combustion_enthalpy():
    eng = ReactionEngine()
    out = eng.predict([Molecule.parse("methane"), Molecule.from_smiles("O=O")],
                      Conditions(temperature_c=600))
    comb = [o for o in out if o.rule_id == "combustion"][0]
    th = comb.thermodynamics()
    # bond-enthalpy estimate for CH4 combustion ≈ -802 kJ/mol
    assert -850 < th["delta_h_kj"] < -750
    assert "exothermic" in th["type_fa"]
    assert "خودبه‌خودی" in th["spontaneity_fa"]


def test_kinetics_yield_and_rate():
    eng = ReactionEngine()
    cond = Conditions(temperature_c=80, catalyst="H2SO4")
    out = eng.predict([Molecule.parse("acetic acid"), Molecule.parse("ethanol")], cond)
    est = [o for o in out if o.rule_id == "esterification"][0]
    kin = est.kinetics()
    assert kin["feasible"] and 0 < kin["estimated_yield_pct"] <= 100
    assert kin["relative_rate"] > 1  # faster than at 25°C


def test_kinetics_infeasible_zero_yield():
    eng = ReactionEngine()
    out = eng.predict([Molecule.parse("acetic acid"), Molecule.parse("ethanol")],
                      Conditions())  # no catalyst/heat
    est = [o for o in out if o.rule_id == "esterification"][0]
    assert est.kinetics()["estimated_yield_pct"] == 0


def test_stoichiometry_limiting_reagent():
    from chemlab.stoichiometry import stoichiometry
    res = stoichiometry(
        [Molecule.parse("methane"), Molecule.from_smiles("O=O")],
        [Molecule.from_smiles("O=C=O"), Molecule.from_smiles("O")],
        [1, 2], [1, 2],
        [{"value": 16, "unit": "g"}, {"value": 64, "unit": "g"}])
    assert res["limiting_reagent"] == "CH4"
    # 1 mol CH4 -> 1 mol CO2 ≈ 44 g
    assert 43 < res["products"][0]["theoretical_g"] < 45


def test_solution_ph():
    from chemlab.solution import ph_of_solution
    assert abs(ph_of_solution(Molecule.parse("hydrochloric acid"), 0.1)["ph"] - 1.0) < 0.1
    assert abs(ph_of_solution(Molecule.parse("acetic acid"), 0.1)["ph"] - 2.87) < 0.2
    assert ph_of_solution(Molecule.parse("sodium hydroxide"), 0.01)["ph"] > 11


def test_distillation_order():
    from chemlab.separation import distill
    res = distill([Molecule.parse("water"), Molecule.parse("ethanol"),
                   Molecule.parse("acetone")])
    order = [f["formula"] for f in res["fractions"]]
    assert order[0] == "C3H6O" and order[-1] == "H2O"  # acetone first, water last


def test_oxidation_states():
    from chemlab.redox import oxidation_states, cell_potential, electrolysis
    co2 = oxidation_states(Molecule.from_smiles("O=C=O"))["by_element"]
    assert co2["C"] == [4] and co2["O"] == [-2]
    cell = cell_potential("Cu2+/Cu", "Zn2+/Zn")
    assert abs(cell["emf"] - 1.10) < 0.01 and cell["spontaneous"]
    assert electrolysis("NaCl")["ok"]


def test_admet_and_interaction():
    from chemlab.admet import admet, interaction
    a = admet(Molecule.parse("caffeine"))
    assert "absorption_fa" in a and "toxicity_fa" in a
    ix = interaction("aspirin", "warfarin")
    assert ix["has_interaction"]


def test_spectra_prediction():
    from chemlab.spectra import full_spectra
    s = full_spectra(Molecule.parse("paracetamol"))
    assert s["nmr_h"]["num_signals"] >= 4
    assert any("آمید" in b["assignment_fa"] for b in s["ir"]["bands"])
    assert abs(s["ms"]["molecular_ion_mz"] - 151.063) < 0.01


def test_nlp_parsing():
    from chemlab.nlp import parse_scenario
    p = parse_scenario("استیک اسید و اتانول با کاتالیزور سولفوریک اسید در ۸۰ درجه")
    assert "اتانول" in p["reactants"]
    assert p["_conditions_obj"].temperature_c == 80
    assert p["_conditions_obj"].catalyst == "H2SO4"


def test_react_nl_end_to_end():
    sp = ScenarioProcessor()
    res = sp.react_nl("متان را در حضور اکسیژن در دمای ۶۰۰ درجه بسوزان")
    assert res["ok"] and res["main_product"]
    assert "→" in res["main_product"]["balanced_equation"]


def test_expanded_reaction_count():
    from chemlab.reaction_data import REACTION_RULES
    assert len(REACTION_RULES) >= 15


def test_qsar_model_trains_and_predicts():
    from chemlab.qsar import demo_model
    m = demo_model()
    assert m.n_train >= 40
    assert m.metrics.get("cv_auc", 0) > 0.6  # real learned signal
    # caffeine (BBB+) should score higher than glucose (BBB-)
    caf = m.predict("Cn1cnc2c1c(=O)n(C)c(=O)n2C")["probability_active"]
    glu = m.predict("OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")["probability_active"]
    assert caf > glu


def test_generator_makes_novel_valid_molecules():
    from chemlab.generator import mutate
    from rdkit import Chem
    seen = set()
    for _ in range(20):
        child = mutate("c1ccccc1")
        if child:
            assert Chem.MolFromSmiles(child) is not None
            seen.add(child)
    assert len(seen) >= 3  # produced several distinct valid novel molecules


def test_discovery_ranks_candidates():
    from chemlab.discovery import discover
    res = discover(objective="qed", population_size=20, generations=3, top_k=5)
    assert res["ok"] and res["candidates"]
    cands = res["candidates"]
    # composite scores are present and sorted descending
    scores = [c["composite_score"] for c in cands]
    assert scores == sorted(scores, reverse=True)
    assert cands[0]["rank"] == 1
    for c in cands:
        assert "activity_score" in c and "admet_score" in c and "synthesizability_score" in c


def test_docking_runs_or_estimates():
    from chemlab.docking import dock
    res = dock("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    assert res["ok"]
    assert "binding_affinity_kcal_mol" in res
    assert res["binding_affinity_kcal_mol"] < 0  # affinity is negative kcal/mol


def test_discovery_with_docking():
    from chemlab.discovery import discover
    res = discover(objective="qed", population_size=14, generations=2,
                   top_k=3, use_docking=True)
    assert res["ok"] and res["candidates"]
    assert all("docking_score" in c for c in res["candidates"])


def test_external_db_graceful_when_blocked():
    from chemlab.external_db import pubchem_by_name
    res = pubchem_by_name("aspirin")
    # network is blocked here; must fail gracefully, never raise
    assert "ok" in res
    if not res["ok"]:
        assert "error_fa" in res


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
