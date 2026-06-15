"""
molecule.py — Molecule model built on top of RDKit.

A `Molecule` wraps an RDKit Mol and exposes a clean, JSON-friendly view of a
chemical species: structure (SMILES / formula), bonds, rings, functional
descriptors, drug-likeness, and a 2D depiction (SVG). Input can be SMILES,
InChI, molecular formula (best-effort), or a known common/Persian name.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors, Crippen, Lipinski
from rdkit.Chem.Draw import rdMolDraw2D

from .known_compounds import resolve_name


BOND_TYPE_FA = {
    Chem.BondType.SINGLE: "یگانه (کووالانسی)",
    Chem.BondType.DOUBLE: "دوگانه",
    Chem.BondType.TRIPLE: "سه‌گانه",
    Chem.BondType.AROMATIC: "آروماتیک",
}


class MoleculeError(ValueError):
    """Raised when an input cannot be parsed into a molecule."""


@dataclass
class Molecule:
    mol: Chem.Mol
    source: str = ""

    # ----- constructors -------------------------------------------------
    @classmethod
    def from_smiles(cls, smiles: str) -> "Molecule":
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise MoleculeError(f"SMILES نامعتبر است: {smiles!r}")
        return cls(mol, source=smiles)

    @classmethod
    def from_inchi(cls, inchi: str) -> "Molecule":
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            raise MoleculeError(f"InChI نامعتبر است: {inchi!r}")
        return cls(mol, source=inchi)

    @classmethod
    def parse(cls, text: str) -> "Molecule":
        """Best-effort parse: name -> SMILES -> InChI."""
        text = text.strip()
        if not text:
            raise MoleculeError("ورودی خالی است")
        # 1) Known common / Persian name.
        smiles = resolve_name(text)
        if smiles:
            return cls.from_smiles(smiles)
        # 2) InChI.
        if text.lower().startswith("inchi="):
            return cls.from_inchi(text)
        # 3) Raw SMILES.
        mol = Chem.MolFromSmiles(text)
        if mol is not None:
            return cls(mol, source=text)
        raise MoleculeError(
            f"نتوانستم {text!r} را تفسیر کنم. SMILES، InChI یا نام شناخته‌شده بدهید."
        )

    # ----- identity -----------------------------------------------------
    @property
    def smiles(self) -> str:
        return Chem.MolToSmiles(self.mol)

    @property
    def formula(self) -> str:
        return rdMolDescriptors.CalcMolFormula(self.mol)

    @property
    def molar_mass(self) -> float:
        return Descriptors.MolWt(self.mol)

    @property
    def inchi(self) -> str:
        return Chem.MolToInchi(self.mol)

    @property
    def inchikey(self) -> str:
        return Chem.InchiToInchiKey(self.inchi)

    # ----- structure ----------------------------------------------------
    def atom_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        mol = Chem.AddHs(self.mol)
        for atom in mol.GetAtoms():
            counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
        return dict(sorted(counts.items()))

    def bonds(self) -> list[dict]:
        out = []
        for b in self.mol.GetBonds():
            out.append({
                "begin": b.GetBeginAtom().GetSymbol(),
                "begin_idx": b.GetBeginAtomIdx(),
                "end": b.GetEndAtom().GetSymbol(),
                "end_idx": b.GetEndAtomIdx(),
                "type": str(b.GetBondType()),
                "type_fa": BOND_TYPE_FA.get(b.GetBondType(), str(b.GetBondType())),
                "is_in_ring": b.IsInRing(),
                "is_conjugated": b.GetIsConjugated(),
            })
        return out

    def rings(self) -> dict:
        ri = self.mol.GetRingInfo()
        return {
            "num_rings": ri.NumRings(),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(self.mol),
            "ring_sizes": [len(r) for r in ri.AtomRings()],
        }

    # ----- physicochemical descriptors ----------------------------------
    def descriptors(self) -> dict:
        m = self.mol
        return {
            "molar_mass": round(Descriptors.MolWt(m), 3),
            "exact_mass": round(Descriptors.ExactMolWt(m), 4),
            "logp": round(Crippen.MolLogP(m), 3),          # lipophilicity
            "tpsa": round(rdMolDescriptors.CalcTPSA(m), 2),  # polar surface area
            "h_bond_donors": Lipinski.NumHDonors(m),
            "h_bond_acceptors": Lipinski.NumHAcceptors(m),
            "rotatable_bonds": Lipinski.NumRotatableBonds(m),
            "heavy_atoms": m.GetNumHeavyAtoms(),
            "formal_charge": Chem.GetFormalCharge(m),
            "num_rings": rdMolDescriptors.CalcNumRings(m),
            "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
        }

    # ----- depiction ----------------------------------------------------
    def to_svg(self, width: int = 360, height: int = 280) -> str:
        mol = Chem.Mol(self.mol)
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:
            pass
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    # ----- export -------------------------------------------------------
    def to_dict(self, include_svg: bool = True) -> dict:
        data = {
            "smiles": self.smiles,
            "formula": self.formula,
            "molar_mass": round(self.molar_mass, 3),
            "inchikey": self.inchikey,
            "atom_counts": self.atom_counts(),
            "bonds": self.bonds(),
            "rings": self.rings(),
            "descriptors": self.descriptors(),
        }
        if include_svg:
            data["svg"] = self.to_svg()
        return data
