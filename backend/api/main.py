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
)
from chemlab.conditions import Conditions, TECHNIQUE_FA, Technique
from chemlab.scenario import _reverse_lookup

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


@app.post("/api/hypothesize")
def hypothesize(req: HypothesisRequest) -> dict:
    return processor.hypothesize(
        req.reactants, req.conditions.to_conditions(), req.pharma_mode, req.top_k
    )


# ----------------------------- static UI --------------------------------
@app.get("/")
def index() -> FileResponse:
    return FileResponse(os.path.join(FRONTEND_DIR, "index.html"))


if os.path.isdir(FRONTEND_DIR):
    app.mount("/static", StaticFiles(directory=FRONTEND_DIR), name="static")
