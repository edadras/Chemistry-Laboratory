"""
generator.py — De novo molecular generation (genetic algorithm).

A graph-based evolutionary generator (in the spirit of GB-GA, Jensen 2019):
starting from seed molecules it applies valence-respecting graph mutations
(add/replace/delete atom, append fragment) and crossover, keeps chemically
valid offspring, scores them with a user-supplied fitness function, and evolves
the population toward higher fitness. No GPU or training needed — it really
produces novel, valid molecules optimised for an objective.

For deep generative models (VAE / GNN / diffusion) you would swap this module
for a trained network; the rest of the discovery pipeline is unchanged.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Callable

from rdkit import Chem
from rdkit.Chem import RWMol
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

# small building blocks appended during mutation
_FRAGMENTS = ["C", "N", "O", "F", "Cl", "CC", "CO", "CN", "C=O", "C#N",
              "C(=O)O", "C(=O)N", "OC", "S", "c1ccccc1", "C1CC1", "n1ccccc1"]
_ELEMENTS = [6, 7, 8, 9, 17, 16]  # C N O F Cl S
_DEFAULT_SEEDS = ["c1ccccc1", "CCO", "CC(=O)O", "C1CCCCC1", "c1ccncc1",
                  "CC(=O)N", "CCN", "O=C(O)c1ccccc1"]


def _free_atoms(mol: Chem.Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetTotalNumHs() > 0]


def _sanitized(mol: Chem.Mol) -> Chem.Mol | None:
    try:
        m = mol.GetMol() if isinstance(mol, RWMol) else mol
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return None


def _mut_add_atom(mol):
    free = _free_atoms(mol)
    if not free:
        return None
    rw = RWMol(mol)
    new = rw.AddAtom(Chem.Atom(random.choice(_ELEMENTS)))
    rw.AddBond(random.choice(free), new, Chem.BondType.SINGLE)
    return _sanitized(rw)


def _mut_replace_atom(mol):
    cand = [a.GetIdx() for a in mol.GetAtoms()
            if not a.GetIsAromatic() and a.GetSymbol() in ("C", "N", "O")]
    if not cand:
        return None
    rw = RWMol(mol)
    rw.GetAtomWithIdx(random.choice(cand)).SetAtomicNum(random.choice([6, 7, 8]))
    return _sanitized(rw)


def _mut_delete_atom(mol):
    terminal = [a.GetIdx() for a in mol.GetAtoms()
                if a.GetDegree() == 1 and not a.GetIsAromatic()]
    if not terminal or mol.GetNumAtoms() <= 3:
        return None
    rw = RWMol(mol)
    rw.RemoveAtom(random.choice(terminal))
    return _sanitized(rw)


def _mut_append_fragment(mol):
    free = _free_atoms(mol)
    if not free:
        return None
    frag = Chem.MolFromSmiles(random.choice(_FRAGMENTS))
    if frag is None:
        return None
    combo = Chem.CombineMols(mol, frag)
    rw = RWMol(combo)
    n1 = mol.GetNumAtoms()
    frag_free = [n1 + a.GetIdx() for a in frag.GetAtoms() if a.GetTotalNumHs() > 0]
    if not frag_free:
        return None
    rw.AddBond(random.choice(free), random.choice(frag_free), Chem.BondType.SINGLE)
    return _sanitized(rw)


_MUTATIONS = [_mut_append_fragment, _mut_add_atom, _mut_replace_atom,
              _mut_delete_atom, _mut_append_fragment]


def mutate(smiles: str, tries: int = 12, max_atoms: int = 40) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for _ in range(tries):
        out = random.choice(_MUTATIONS)(mol)
        if out is None or out.GetNumAtoms() > max_atoms:
            continue
        smi = Chem.MolToSmiles(out)
        if smi and smi != smiles:
            return smi
    return None


def crossover(smiles_a: str, smiles_b: str) -> str | None:
    """Swap a terminal fragment of A with one from B (acyclic single-bond cut)."""
    a, b = Chem.MolFromSmiles(smiles_a), Chem.MolFromSmiles(smiles_b)
    if a is None or b is None:
        return None
    fa, fb = _free_atoms(a), _free_atoms(b)
    if not fa or not fb:
        return None
    combo = Chem.CombineMols(a, b)
    rw = RWMol(combo)
    rw.AddBond(random.choice(fa), a.GetNumAtoms() + random.choice(fb),
               Chem.BondType.SINGLE)
    out = _sanitized(rw)
    return Chem.MolToSmiles(out) if out is not None else None


@dataclass
class GAResult:
    population: list[tuple[str, float]]  # (smiles, fitness) sorted desc
    generations: int
    n_evaluated: int


def evolve(fitness: Callable[[str], float], seeds: list[str] | None = None,
           population_size: int = 40, generations: int = 8,
           elite: int = 8, seed: int | None = 42) -> GAResult:
    """Run the genetic algorithm and return the ranked final population."""
    if seed is not None:
        random.seed(seed)
    seeds = seeds or _DEFAULT_SEEDS
    pop = {s: fitness(s) for s in seeds if Chem.MolFromSmiles(s)}
    evaluated = set(pop)

    for _ in range(generations):
        ranked = sorted(pop.items(), key=lambda kv: kv[1], reverse=True)
        parents = [s for s, _ in ranked[:max(elite, 2)]]
        children: dict[str, float] = {}
        while len(children) < population_size:
            if len(parents) >= 2 and random.random() < 0.35:
                child = crossover(random.choice(parents), random.choice(parents))
            else:
                child = mutate(random.choice(parents))
            if not child or child in evaluated or child in children:
                if len(evaluated) > 5000:
                    break
                continue
            children[child] = fitness(child)
            evaluated.add(child)
        # next generation = elites + children
        merged = dict(ranked[:elite])
        merged.update(children)
        pop = merged

    ranked = sorted(pop.items(), key=lambda kv: kv[1], reverse=True)
    return GAResult(population=ranked, generations=generations,
                    n_evaluated=len(evaluated))
