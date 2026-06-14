"""
known_compounds.py — Dictionary of well-known compounds and drugs.

Maps human-friendly names (English + Persian, plus a few synonyms) to
canonical SMILES so users can type "آسپرین" or "water" instead of a SMILES
string. Also powers the pharmaceutical knowledge base.
"""
from __future__ import annotations

# name (lowercase) -> SMILES
# Each entry registered under English name, Persian name and common synonyms.
_COMPOUNDS: dict[str, dict] = {
    # --- common inorganic / small molecules ---
    "water": {"smiles": "O", "fa": "آب", "syn": ["h2o", "aqua"]},
    "hydrogen peroxide": {"smiles": "OO", "fa": "آب اکسیژنه", "syn": ["h2o2", "peroxide"]},
    "carbon dioxide": {"smiles": "O=C=O", "fa": "دی‌اکسید کربن", "syn": ["co2"]},
    "carbon monoxide": {"smiles": "[C-]#[O+]", "fa": "مونوکسید کربن", "syn": ["co"]},
    "ammonia": {"smiles": "N", "fa": "آمونیاک", "syn": ["nh3"]},
    "methane": {"smiles": "C", "fa": "متان", "syn": ["ch4"]},
    "ethane": {"smiles": "CC", "fa": "اتان", "syn": []},
    "propane": {"smiles": "CCC", "fa": "پروپان", "syn": []},
    "butane": {"smiles": "CCCC", "fa": "بوتان", "syn": []},
    "ethylene": {"smiles": "C=C", "fa": "اتیلن", "syn": ["ethene"]},
    "acetylene": {"smiles": "C#C", "fa": "استیلن", "syn": ["ethyne"]},
    "benzene": {"smiles": "c1ccccc1", "fa": "بنزن", "syn": []},
    "toluene": {"smiles": "Cc1ccccc1", "fa": "تولوئن", "syn": []},
    "phenol": {"smiles": "Oc1ccccc1", "fa": "فنل", "syn": []},
    "methanol": {"smiles": "CO", "fa": "متانول", "syn": ["methyl alcohol"]},
    "ethanol": {"smiles": "CCO", "fa": "اتانول", "syn": ["alcohol", "الکل", "ethyl alcohol"]},
    "isopropanol": {"smiles": "CC(C)O", "fa": "ایزوپروپانول", "syn": ["ipa", "rubbing alcohol"]},
    "glycerol": {"smiles": "OCC(O)CO", "fa": "گلیسرول", "syn": ["glycerin", "گلیسیرین"]},
    "acetic acid": {"smiles": "CC(=O)O", "fa": "اسید استیک", "syn": ["ethanoic acid", "vinegar"]},
    "formic acid": {"smiles": "C(=O)O", "fa": "اسید فرمیک", "syn": ["methanoic acid"]},
    "acetone": {"smiles": "CC(=O)C", "fa": "استون", "syn": ["propanone"]},
    "formaldehyde": {"smiles": "C=O", "fa": "فرمالدهید", "syn": ["methanal"]},
    "acetaldehyde": {"smiles": "CC=O", "fa": "استالدهید", "syn": ["ethanal"]},
    "glucose": {"smiles": "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "fa": "گلوکز", "syn": ["dextrose"]},
    "sucrose": {"smiles": "OC[C@H]1O[C@@](CO)(O[C@@H]2OC(CO)[C@@H](O)[C@H]2O)[C@H](O)[C@@H]1O",
                "fa": "ساکارز", "syn": ["sugar", "شکر", "قند"]},
    "sodium chloride": {"smiles": "[Na+].[Cl-]", "fa": "کلرید سدیم", "syn": ["nacl", "salt", "نمک"]},
    "sodium hydroxide": {"smiles": "[Na+].[OH-]", "fa": "سود سوزآور", "syn": ["naoh", "lye", "caustic soda"]},
    "sulfuric acid": {"smiles": "OS(=O)(=O)O", "fa": "اسید سولفوریک", "syn": ["h2so4"]},
    "nitric acid": {"smiles": "O[N+](=O)[O-]", "fa": "اسید نیتریک", "syn": ["hno3"]},
    "hydrochloric acid": {"smiles": "Cl", "fa": "اسید کلریدریک", "syn": ["hcl", "جوهر نمک"]},
    "sodium bicarbonate": {"smiles": "C(=O)(O)[O-].[Na+]", "fa": "جوش‌شیرین",
                           "syn": ["baking soda", "nahco3", "بی‌کربنات سدیم"]},
    "urea": {"smiles": "NC(=O)N", "fa": "اوره", "syn": []},

    # --- drugs / pharma ---
    "aspirin": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "fa": "آسپرین",
                "syn": ["acetylsalicylic acid", "استیل سالیسیلیک اسید"]},
    "salicylic acid": {"smiles": "O=C(O)c1ccccc1O", "fa": "سالیسیلیک اسید", "syn": []},
    "paracetamol": {"smiles": "CC(=O)Nc1ccc(O)cc1", "fa": "استامینوفن",
                    "syn": ["acetaminophen", "tylenol", "استامینوفن"]},
    "ibuprofen": {"smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "fa": "ایبوپروفن", "syn": ["advil"]},
    "caffeine": {"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "fa": "کافئین", "syn": []},
    "penicillin g": {"smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
                     "fa": "پنی‌سیلین جی", "syn": ["benzylpenicillin", "penicillin", "پنی‌سیلین"]},
    "amoxicillin": {"smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O",
                    "fa": "آموکسی‌سیلین", "syn": ["amoxil"]},
    "ascorbic acid": {"smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "fa": "ویتامین ث",
                      "syn": ["vitamin c", "ویتامین سی", "vit c"]},
    "morphine": {"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
                 "fa": "مورفین", "syn": []},
    "nicotine": {"smiles": "CN1CCC[C@H]1c1cccnc1", "fa": "نیکوتین", "syn": []},
    "dopamine": {"smiles": "NCCc1ccc(O)c(O)c1", "fa": "دوپامین", "syn": []},
    "adrenaline": {"smiles": "CNC[C@H](O)c1ccc(O)c(O)c1", "fa": "آدرنالین", "syn": ["epinephrine"]},
    "penicillin v": {"smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)COc3ccccc3)C(=O)N2[C@H]1C(=O)O",
                     "fa": "پنی‌سیلین وی", "syn": ["phenoxymethylpenicillin"]},
    "metformin": {"smiles": "CN(C)C(=N)N=C(N)N", "fa": "متفورمین", "syn": []},
    "diazepam": {"smiles": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", "fa": "دیازپام", "syn": ["valium"]},
}


def _build_index() -> dict[str, str]:
    idx: dict[str, str] = {}
    for canonical, info in _COMPOUNDS.items():
        smiles = info["smiles"]
        idx[canonical.lower()] = smiles
        if info.get("fa"):
            idx[info["fa"].strip().lower()] = smiles
        for syn in info.get("syn", []):
            idx[syn.strip().lower()] = smiles
    return idx


_INDEX = _build_index()


def resolve_name(text: str) -> str | None:
    """Return SMILES for a known compound name, else None."""
    return _INDEX.get(text.strip().lower())


def all_compounds() -> list[dict]:
    """List the catalogue (for the UI)."""
    return [
        {"name": name, "name_fa": info.get("fa"), "smiles": info["smiles"]}
        for name, info in _COMPOUNDS.items()
    ]
