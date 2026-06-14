"""
admet.py — Deeper pharmacology: ADMET-style predictions & drug interactions.

Rule-of-thumb estimates (transparent, descriptor-based) for:
  * Absorption — GI absorption & oral bioavailability (Veber / Egan).
  * Distribution — blood–brain-barrier permeability.
  * Metabolism / lipophilicity flags.
  * Excretion — rough renal vs hepatic tendency.
  * Toxicity — structural-alert toxicophores.
Plus a small known drug–drug interaction table. These are educational
estimates, NOT clinical advice.
"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

from .molecule import Molecule
from .pharma import STRUCTURAL_ALERTS
from .known_compounds import drug_info


def admet(mol: Molecule) -> dict:
    m = mol.mol
    mw = Descriptors.MolWt(m)
    logp = Crippen.MolLogP(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)
    rotb = Lipinski.NumRotatableBonds(m)

    # Absorption
    gi_high = tpsa <= 140 and rotb <= 10 and mw <= 500
    absorption = "جذب گوارشی بالا" if gi_high else "جذب گوارشی پایین"

    # Distribution — blood-brain barrier
    bbb = (tpsa <= 90 and mw <= 450 and 1 <= logp <= 4)
    distribution = "احتمال عبور از سد خونی-مغزی (اثر مرکزی)" if bbb \
        else "بعید است از سد خونی-مغزی عبور کند (اثر محیطی)"

    # Metabolism — lipophilicity tends to raise CYP metabolism
    if logp > 3:
        metabolism = "چربی‌دوست؛ احتمال متابولیسم کبدی (CYP) بالا"
    elif logp < 0:
        metabolism = "آب‌دوست؛ متابولیسم کم، دفع کلیوی محتمل"
    else:
        metabolism = "متابولیسم متعادل"

    # Excretion — small polar molecules → renal
    excretion = "دفع کلیوی (مولکول کوچک/قطبی)" if (mw < 350 and logp < 1) \
        else "دفع کبدی/صفراوی محتمل"

    # Toxicity — structural alerts
    alerts = []
    for smarts, label in STRUCTURAL_ALERTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and m.HasSubstructMatch(patt):
            alerts.append(label)
    tox = "هشدار ساختاری دارد" if alerts else "بدون هشدار ساختاری شناخته‌شده"

    return {
        "absorption_fa": absorption,
        "gi_absorption_high": gi_high,
        "distribution_fa": distribution,
        "bbb_permeant": bbb,
        "metabolism_fa": metabolism,
        "excretion_fa": excretion,
        "toxicity_fa": tox,
        "toxicophores": alerts,
        "descriptors": {"mw": round(mw, 1), "logp": round(logp, 2),
                        "tpsa": round(tpsa, 1), "hbd": hbd, "hba": hba, "rotb": rotb},
        "is_estimate": True,
        "disclaimer_fa": "تخمین آموزشی بر پایه‌ی توصیف‌گرها — جایگزین داده‌ی بالینی نیست.",
    }


# Known drug–drug interactions (canonical English names).
_INTERACTIONS = {
    frozenset({"aspirin", "warfarin"}):
        "افزایش خطر خونریزی (هر دو روی انعقاد اثر می‌گذارند)",
    frozenset({"aspirin", "ibuprofen"}):
        "رقابت بر سر COX؛ ایبوپروفن می‌تواند اثر ضدپلاکتی آسپرین را کاهش دهد",
    frozenset({"warfarin", "diclofenac"}):
        "افزایش خطر خونریزی گوارشی",
    frozenset({"fluoxetine", "tramadol"}):
        "خطر سندرم سروتونین",
    frozenset({"ciprofloxacin", "metformin"}):
        "احتمال تغییر قند خون؛ پایش لازم است",
    frozenset({"omeprazole", "diazepam"}):
        "امپرازول متابولیسم دیازپام را کند می‌کند (افزایش اثر)",
}


def interaction(name_a: str, name_b: str) -> dict:
    a = drug_info(name_a)
    b = drug_info(name_b)
    na = a["name"] if a else name_a.strip().lower()
    nb = b["name"] if b else name_b.strip().lower()
    note = _INTERACTIONS.get(frozenset({na, nb}))
    return {
        "drug_a": na, "drug_b": nb,
        "interaction_fa": note or "تداخل شناخته‌شده‌ای در پایگاه ثبت نشده است",
        "has_interaction": note is not None,
    }
