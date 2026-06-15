"""
FastAPI backend for the Chemistry Laboratory.

Exposes the `chemlab` engine over HTTP and serves the single-page web UI.

Run:  uvicorn api.main:app --reload   (from the backend/ directory)
"""
from __future__ import annotations

import os
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from chemlab import (
    ScenarioProcessor, Molecule, MoleculeError, periodic_table,
    all_compounds, DrugProfile, indication_for,
    ph_of_solution, solubility, distill, oxidation_states, cell_potential,
    electrolysis, STANDARD_POTENTIALS, admet, interaction, full_spectra,
)
from chemlab.conditions import Conditions, TECHNIQUE_FA, Technique
from chemlab.scenario import _reverse_lookup
from chemlab.known_compounds import all_drugs
from chemlab.qsar import demo_model
from chemlab.discovery import discover
from chemlab.docking import dock as run_dock
from chemlab import external_db


def _parse_or_404(query: str) -> Molecule:
    try:
        return Molecule.parse(query)
    except MoleculeError as e:
        raise HTTPException(status_code=400, detail=str(e))

app = FastAPI(title="آزمایشگاه شیمی هوشمند", version="0.1.0")
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"],
)

processor = ScenarioProcessor()

FRONTEND_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "frontend")


# ----------------------------- schemas ----------------------------------
class ConditionsIn(BaseModel):
    temperature_c: float = 25.0
    pressure_atm: float = 1.0
    humidity_pct: float = 50.0
    ph: Optional[float] = None
    catalyst: Optional[str] = None
    solvent: Optional[str] = None
    light: bool = False
    duration_min: Optional[float] = None
    techniques: list[str] = Field(default_factory=list)

    def to_conditions(self) -> Conditions:
        return Conditions.from_dict(self.model_dump())


class ReactRequest(BaseModel):
    reactants: list[str]
    conditions: ConditionsIn = Field(default_factory=ConditionsIn)
    pharma_mode: bool = False


class ReverseRequest(BaseModel):
    target: str
    max_depth: int = 3


class SynthesizeRequest(BaseModel):
    target: str
    multistep: bool = False
    max_steps: int = 4


class NLRequest(BaseModel):
    text: str
    pharma_mode: bool = False


class QuantRequest(BaseModel):
    reactants: list[str]
    conditions: ConditionsIn = Field(default_factory=ConditionsIn)
    amounts: dict[str, dict] = Field(default_factory=dict)
    actual_yield_g: Optional[float] = None


class PHRequest(BaseModel):
    compound: str
    concentration_m: float = 0.1


class DistillRequest(BaseModel):
    components: list[str]


class CellRequest(BaseModel):
    cathode: str
    anode: str


class InteractionRequest(BaseModel):
    drug_a: str
    drug_b: str


class DiscoverRequest(BaseModel):
    objective: str = "activity"           # "activity" or "qed"
    seeds: Optional[list[str]] = None
    population_size: int = 30
    generations: int = 6
    top_k: int = 10
    use_docking: bool = False


class QSARPredictRequest(BaseModel):
    smiles: str


class DockRequest(BaseModel):
    smiles: str
    exhaustiveness: int = 4


class HypothesisRequest(BaseModel):
    reactants: list[str]
    conditions: ConditionsIn = Field(default_factory=ConditionsIn)
    pharma_mode: bool = True
    top_k: int = 10


# ----------------------------- endpoints --------------------------------
@app.get("/api/health")
def health() -> dict:
    return {"ok": True, "service": "chemistry-laboratory", "version": app.version}


@app.get("/api/techniques")
def techniques() -> dict:
    return {"techniques": [{"id": t.value, "name_fa": TECHNIQUE_FA[t]} for t in Technique]}


@app.get("/api/compounds")
def compounds() -> dict:
    return {"compounds": all_compounds()}


@app.get("/api/drugs")
def drugs() -> dict:
    return {"drugs": all_drugs()}


@app.get("/api/elements")
def elements() -> dict:
    return {"elements": [e.to_dict() for e in periodic_table.all()]}


@app.get("/api/element/{query}")
def element(query: str) -> dict:
    try:
        if query.isdigit():
            el = periodic_table.by_number(int(query))
        else:
            try:
                el = periodic_table.by_symbol(query)
            except ValueError:
                el = periodic_table.by_name(query)
        return el.to_dict()
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))


@app.get("/api/molecule/{query:path}")
def molecule(query: str) -> dict:
    try:
        mol = Molecule.parse(query)
    except MoleculeError as e:
        raise HTTPException(status_code=400, detail=str(e))
    data = mol.to_dict()
    name = _reverse_lookup(mol.smiles)
    data["pharma"] = DrugProfile(mol).to_dict()
    if name:
        data["known_name"] = name
        data["indication"] = indication_for(name)
    return data


@app.get("/api/molecule_svg/{query:path}")
def molecule_svg(query: str) -> Response:
    try:
        mol = Molecule.parse(query)
    except MoleculeError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return Response(content=mol.to_svg(), media_type="image/svg+xml")


@app.post("/api/react")
def react(req: ReactRequest) -> dict:
    return processor.react(req.reactants, req.conditions.to_conditions(), req.pharma_mode)


@app.post("/api/reverse")
def reverse(req: ReverseRequest) -> dict:
    return processor.reverse(req.target, req.max_depth)


@app.post("/api/synthesize")
def synthesize(req: SynthesizeRequest) -> dict:
    return processor.synthesize(req.target, req.multistep, req.max_steps)


@app.post("/api/hypothesize")
def hypothesize(req: HypothesisRequest) -> dict:
    return processor.hypothesize(
        req.reactants, req.conditions.to_conditions(), req.pharma_mode, req.top_k
    )


@app.post("/api/scenario_nl")
def scenario_nl(req: NLRequest) -> dict:
    return processor.react_nl(req.text, req.pharma_mode)


@app.post("/api/quantitative")
def quantitative(req: QuantRequest) -> dict:
    return processor.react_quantitative(
        req.reactants, req.conditions.to_conditions(), req.amounts, req.actual_yield_g)


@app.post("/api/ph")
def ph(req: PHRequest) -> dict:
    return ph_of_solution(_parse_or_404(req.compound), req.concentration_m)


@app.get("/api/solubility/{query:path}")
def solubility_ep(query: str) -> dict:
    return solubility(_parse_or_404(query))


@app.post("/api/distill")
def distill_ep(req: DistillRequest) -> dict:
    mols = [_parse_or_404(c) for c in req.components if c.strip()]
    return distill(mols)


@app.get("/api/oxidation/{query:path}")
def oxidation_ep(query: str) -> dict:
    return oxidation_states(_parse_or_404(query))


@app.post("/api/cell")
def cell_ep(req: CellRequest) -> dict:
    return cell_potential(req.cathode, req.anode)


@app.get("/api/half_reactions")
def half_reactions() -> dict:
    return {"half_reactions": [{"id": k, "e": v["e"], "fa": v["fa"]}
                               for k, v in STANDARD_POTENTIALS.items()]}


@app.get("/api/electrolysis/{electrolyte}")
def electrolysis_ep(electrolyte: str) -> dict:
    return electrolysis(electrolyte)


@app.get("/api/admet/{query:path}")
def admet_ep(query: str) -> dict:
    return admet(_parse_or_404(query))


@app.post("/api/interaction")
def interaction_ep(req: InteractionRequest) -> dict:
    return interaction(req.drug_a, req.drug_b)


@app.get("/api/spectra/{query:path}")
def spectra_ep(query: str) -> dict:
    return full_spectra(_parse_or_404(query))


# ----- drug discovery (ML closed loop) ----------------------------------
@app.get("/api/qsar_model")
def qsar_model() -> dict:
    return demo_model().to_dict()


@app.post("/api/qsar_predict")
def qsar_predict(req: QSARPredictRequest) -> dict:
    return demo_model().predict(req.smiles)


@app.post("/api/discover")
def discover_ep(req: DiscoverRequest) -> dict:
    return discover(objective=req.objective, seeds=req.seeds,
                    population_size=min(req.population_size, 60),
                    generations=min(req.generations, 12), top_k=req.top_k,
                    use_docking=req.use_docking)


@app.post("/api/dock")
def dock_ep(req: DockRequest) -> dict:
    return run_dock(req.smiles, exhaustiveness=min(req.exhaustiveness, 8))


@app.get("/api/pubchem/{name:path}")
def pubchem_ep(name: str) -> dict:
    return external_db.pubchem_by_name(name)


@app.get("/api/chembl/{chembl_id}")
def chembl_ep(chembl_id: str) -> dict:
    return external_db.chembl_molecule(chembl_id)


# ----------------------------- static UI --------------------------------
@app.get("/")
def index() -> FileResponse:
    return FileResponse(os.path.join(FRONTEND_DIR, "index.html"))


if os.path.isdir(FRONTEND_DIR):
    app.mount("/static", StaticFiles(directory=FRONTEND_DIR), name="static")
