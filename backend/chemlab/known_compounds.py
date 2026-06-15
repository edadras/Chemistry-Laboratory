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
    "oxygen": {"smiles": "O=O", "fa": "اکسیژن", "syn": ["o2", "dioxygen"]},
    "hydrogen": {"smiles": "[H][H]", "fa": "هیدروژن", "syn": ["h2", "dihydrogen"]},
    "nitrogen": {"smiles": "N#N", "fa": "نیتروژن", "syn": ["n2", "dinitrogen"]},
    "chlorine": {"smiles": "ClCl", "fa": "کلر", "syn": ["cl2", "گاز کلر"]},
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
    "acetic acid": {"smiles": "CC(=O)O", "fa": "اسید استیک",
                    "syn": ["ethanoic acid", "vinegar", "استیک اسید"]},
    "formic acid": {"smiles": "C(=O)O", "fa": "اسید فرمیک",
                    "syn": ["methanoic acid", "فرمیک اسید"]},
    "acetone": {"smiles": "CC(=O)C", "fa": "استون", "syn": ["propanone"]},
    "formaldehyde": {"smiles": "C=O", "fa": "فرمالدهید", "syn": ["methanal"]},
    "acetaldehyde": {"smiles": "CC=O", "fa": "استالدهید", "syn": ["ethanal"]},
    "glucose": {"smiles": "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "fa": "گلوکز", "syn": ["dextrose"]},
    "sucrose": {"smiles": "OC[C@H]1O[C@@](CO)(O[C@@H]2OC(CO)[C@@H](O)[C@H]2O)[C@H](O)[C@@H]1O",
                "fa": "ساکارز", "syn": ["sugar", "شکر", "قند"]},
    "sodium chloride": {"smiles": "[Na+].[Cl-]", "fa": "کلرید سدیم", "syn": ["nacl", "salt", "نمک"]},
    "sodium hydroxide": {"smiles": "[Na+].[OH-]", "fa": "سود سوزآور", "syn": ["naoh", "lye", "caustic soda"]},
    "sulfuric acid": {"smiles": "OS(=O)(=O)O", "fa": "اسید سولفوریک",
                      "syn": ["h2so4", "سولفوریک اسید"]},
    "nitric acid": {"smiles": "O[N+](=O)[O-]", "fa": "اسید نیتریک",
                    "syn": ["hno3", "نیتریک اسید"]},
    "hydrochloric acid": {"smiles": "Cl", "fa": "اسید کلریدریک",
                          "syn": ["hcl", "جوهر نمک", "کلریدریک اسید", "هیدروکلریک اسید"]},
    "sodium bicarbonate": {"smiles": "C(=O)(O)[O-].[Na+]", "fa": "جوش‌شیرین",
                           "syn": ["baking soda", "nahco3", "بی‌کربنات سدیم"]},
    "urea": {"smiles": "NC(=O)N", "fa": "اوره", "syn": []},
    "salicylic acid": {"smiles": "O=C(O)c1ccccc1O", "fa": "سالیسیلیک اسید", "syn": []},
    "nicotine": {"smiles": "CN1CCC[C@H]1c1cccnc1", "fa": "نیکوتین", "syn": []},
    "dopamine": {"smiles": "NCCc1ccc(O)c(O)c1", "fa": "دوپامین", "syn": []},
    "para-aminophenol": {"smiles": "Nc1ccc(O)cc1", "fa": "پاراآمینوفنول", "syn": ["4-aminophenol"]},
}

# ---------------------------------------------------------------------------
# Drug registry — each entry carries therapeutic metadata so the engine can
# report indication (کاربرد درمانی) and class, and group drugs by category.
# Keys: smiles, fa, syn, ind (کاربرد درمانی), cls (دسته دارویی), cat (گروه درمانی)
# ---------------------------------------------------------------------------
_DRUGS: dict[str, dict] = {
    # ---- مسکن‌ها و ضدالتهاب‌ها (analgesics / NSAIDs) ----
    "aspirin": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "fa": "آسپرین",
                "syn": ["acetylsalicylic acid", "استیل سالیسیلیک اسید"],
                "ind": "مسکن، تب‌بر، ضدالتهاب و رقیق‌کننده خون", "cls": "NSAID", "cat": "analgesic"},
    "paracetamol": {"smiles": "CC(=O)Nc1ccc(O)cc1", "fa": "استامینوفن",
                    "syn": ["acetaminophen", "tylenol"],
                    "ind": "تب‌بر و مسکن", "cls": "ضددرد/تب‌بر", "cat": "analgesic"},
    "ibuprofen": {"smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "fa": "ایبوپروفن", "syn": ["advil"],
                  "ind": "ضدالتهاب غیراستروئیدی و مسکن", "cls": "NSAID", "cat": "analgesic"},
    "naproxen": {"smiles": "COc1ccc2cc([C@@H](C)C(=O)O)ccc2c1", "fa": "ناپروکسن", "syn": [],
                 "ind": "ضدالتهاب و مسکن طولانی‌اثر", "cls": "NSAID", "cat": "analgesic"},
    "diclofenac": {"smiles": "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "fa": "دیکلوفناک", "syn": [],
                   "ind": "ضدالتهاب و مسکن (دردهای عضلانی-مفصلی)", "cls": "NSAID", "cat": "analgesic"},
    "ketoprofen": {"smiles": "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1", "fa": "کتوپروفن", "syn": [],
                   "ind": "ضدالتهاب و مسکن", "cls": "NSAID", "cat": "analgesic"},
    "morphine": {"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
                 "fa": "مورفین", "syn": [],
                 "ind": "مسکن قوی برای درد شدید", "cls": "اوپیوئید", "cat": "analgesic"},
    "tramadol": {"smiles": "CN(C)C[C@@H]1CCCC[C@]1(O)c1cccc(OC)c1", "fa": "ترامادول", "syn": [],
                 "ind": "مسکن متوسط تا شدید", "cls": "اوپیوئید ضعیف", "cat": "analgesic"},

    # ---- آنتی‌بیوتیک‌ها (antibiotics) ----
    "penicillin g": {"smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
                     "fa": "پنی‌سیلین جی", "syn": ["benzylpenicillin", "penicillin", "پنی‌سیلین"],
                     "ind": "عفونت‌های باکتریایی گرم‌مثبت", "cls": "بتالاکتام", "cat": "antibiotic"},
    "amoxicillin": {"smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O",
                    "fa": "آموکسی‌سیلین", "syn": ["amoxil"],
                    "ind": "آنتی‌بیوتیک طیف‌وسیع", "cls": "پنی‌سیلین (بتالاکتام)", "cat": "antibiotic"},
    "ampicillin": {"smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
                   "fa": "آمپی‌سیلین", "syn": [],
                   "ind": "آنتی‌بیوتیک طیف‌وسیع", "cls": "پنی‌سیلین (بتالاکتام)", "cat": "antibiotic"},
    "ciprofloxacin": {"smiles": "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
                      "fa": "سیپروفلوکساسین", "syn": ["cipro"],
                      "ind": "عفونت‌های ادراری و تنفسی", "cls": "فلوروکینولون", "cat": "antibiotic"},
    "trimethoprim": {"smiles": "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", "fa": "تری‌متوپریم", "syn": [],
                     "ind": "عفونت ادراری (مهارکننده فولات)", "cls": "ضدفولات", "cat": "antibiotic"},
    "sulfamethoxazole": {"smiles": "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "fa": "سولفامتوکسازول",
                         "syn": [], "ind": "آنتی‌بیوتیک سولفونامیدی", "cls": "سولفونامید",
                         "cat": "antibiotic"},
    "chloramphenicol": {"smiles": "OC[C@@H](NC(=O)C(Cl)Cl)[C@H](O)c1ccc([N+](=O)[O-])cc1",
                        "fa": "کلرامفنیکل", "syn": [],
                        "ind": "آنتی‌بیوتیک طیف‌وسیع", "cls": "آمفنیکل", "cat": "antibiotic"},
    "metronidazole": {"smiles": "Cc1ncc([N+](=O)[O-])n1CCO", "fa": "مترونیدازول", "syn": [],
                      "ind": "عفونت‌های بی‌هوازی و انگلی", "cls": "نیتروایمیدازول", "cat": "antibiotic"},

    # ---- ضدویروس‌ها (antivirals) ----
    "acyclovir": {"smiles": "Nc1nc2n(COCCO)cnc2c(=O)[nH]1", "fa": "آسیکلوویر", "syn": ["aciclovir"],
                  "ind": "عفونت‌های تبخال (هرپس)", "cls": "آنالوگ نوکلئوزیدی", "cat": "antiviral"},
    "zidovudine": {"smiles": "Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O",
                   "fa": "زیدوودین", "syn": ["azt"],
                   "ind": "درمان HIV/ایدز", "cls": "مهارکننده ترانس‌کریپتاز معکوس", "cat": "antiviral"},
    "oseltamivir": {"smiles": "CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H](N)C1",
                    "fa": "اوسلتامیویر", "syn": ["tamiflu"],
                    "ind": "آنفلوانزا (مهار نورامینیداز)", "cls": "مهارکننده نورامینیداز",
                    "cat": "antiviral"},

    # ---- سایر (other therapeutic) ----
    "metformin": {"smiles": "CN(C)C(=N)N=C(N)N", "fa": "متفورمین", "syn": [],
                  "ind": "کنترل قند خون در دیابت نوع ۲", "cls": "بیگوانید", "cat": "antidiabetic"},
    "diazepam": {"smiles": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", "fa": "دیازپام", "syn": ["valium"],
                 "ind": "ضداضطراب و آرام‌بخش", "cls": "بنزودیازپین", "cat": "cns"},
    "caffeine": {"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "fa": "کافئین", "syn": [],
                 "ind": "محرک سیستم عصبی مرکزی", "cls": "محرک گزانتین", "cat": "cns"},
    "ascorbic acid": {"smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "fa": "ویتامین ث",
                      "syn": ["vitamin c", "ویتامین سی", "vit c"],
                      "ind": "ویتامین C، آنتی‌اکسیدان", "cls": "ویتامین", "cat": "vitamin"},
    "adrenaline": {"smiles": "CNC[C@H](O)c1ccc(O)c(O)c1", "fa": "آدرنالین", "syn": ["epinephrine"],
                   "ind": "شوک آنافیلاکسی و ایست قلبی", "cls": "کاتکول‌آمین", "cat": "cardio"},
    "phenylephrine": {"smiles": "CNC[C@H](O)c1cccc(O)c1", "fa": "فنیل‌افرین", "syn": [],
                      "ind": "ضداحتقان بینی و افزاینده فشار خون", "cls": "آگونیست آلفا-آدرنرژیک",
                      "cat": "decongestant"},
    "pseudoephedrine": {"smiles": "CN[C@@H](C)[C@H](O)c1ccccc1", "fa": "سودوافدرین", "syn": [],
                        "ind": "ضداحتقان بینی", "cls": "سمپاتومیمتیک", "cat": "decongestant"},
    "salbutamol": {"smiles": "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "fa": "سالبوتامول", "syn": ["albuterol"],
                   "ind": "گشادکننده برونش (آسم)", "cls": "آگونیست بتا-۲", "cat": "respiratory"},
    "omeprazole": {"smiles": "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1", "fa": "امپرازول", "syn": [],
                   "ind": "زخم معده و رفلاکس (مهار پمپ پروتون)", "cls": "PPI", "cat": "git"},
    "ranitidine": {"smiles": "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1", "fa": "رانیتیدین",
                   "syn": ["zantac"], "ind": "زخم معده (مسدودکننده H2)", "cls": "آنتاگونیست H2",
                   "cat": "git"},

    # ---- آنتی‌هیستامین‌ها (antihistamines) ----
    "diphenhydramine": {"smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1", "fa": "دیفن‌هیدرامین",
                        "syn": ["benadryl"], "ind": "ضدحساسیت و خواب‌آور", "cls": "آنتی‌هیستامین",
                        "cat": "antihistamine"},
    "loratadine": {"smiles": "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1", "fa": "لوراتادین",
                   "syn": ["claritin"], "ind": "ضدحساسیت بدون خواب‌آوری", "cls": "آنتی‌هیستامین",
                   "cat": "antihistamine"},
    "cetirizine": {"smiles": "OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1", "fa": "ستیریزین",
                   "syn": ["zyrtec"], "ind": "ضدحساسیت", "cls": "آنتی‌هیستامین", "cat": "antihistamine"},

    # ---- اعصاب و روان (CNS / antidepressants) ----
    "fluoxetine": {"smiles": "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "fa": "فلوکستین",
                   "syn": ["prozac"], "ind": "ضدافسردگی (مهار بازجذب سروتونین)", "cls": "SSRI",
                   "cat": "cns"},
    "sertraline": {"smiles": "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21", "fa": "سرترالین",
                   "syn": ["zoloft"], "ind": "ضدافسردگی", "cls": "SSRI", "cat": "cns"},
    "gabapentin": {"smiles": "NCC1(CC(=O)O)CCCCC1", "fa": "گاباپنتین", "syn": [],
                   "ind": "درد عصبی و تشنج", "cls": "ضدتشنج", "cat": "cns"},

    # ---- قلبی-عروقی (cardiovascular) ----
    "warfarin": {"smiles": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "fa": "وارفارین",
                 "syn": ["coumadin"], "ind": "ضدانعقاد خون", "cls": "آنتاگونیست ویتامین K",
                 "cat": "cardio"},
    "losartan": {"smiles": "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1", "fa": "لوزارتان",
                 "syn": [], "ind": "فشار خون بالا (آنتاگونیست آنژیوتانسین)", "cls": "ARB",
                 "cat": "cardio"},
    "amlodipine": {"smiles": "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl", "fa": "آملودیپین",
                   "syn": [], "ind": "فشار خون و آنژین (مسدودکننده کانال کلسیم)",
                   "cls": "CCB", "cat": "cardio"},
    "furosemide": {"smiles": "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl", "fa": "فوروزماید",
                   "syn": ["lasix"], "ind": "ادرارآور قوی (نارسایی قلبی/ادم)", "cls": "دیورتیک لوپ",
                   "cat": "cardio"},
    "hydrochlorothiazide": {"smiles": "NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O", "fa": "هیدروکلروتیازید",
                            "syn": ["hctz"], "ind": "ادرارآور و کاهش فشار خون", "cls": "دیورتیک تیازیدی",
                            "cat": "cardio"},

    # ---- ضدمالاریا (antimalarials) ----
    "chloroquine": {"smiles": "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "fa": "کلروکین", "syn": [],
                    "ind": "درمان و پیشگیری مالاریا", "cls": "آمینوکینولین", "cat": "antimalarial"},
    "hydroxychloroquine": {"smiles": "CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "fa": "هیدروکسی‌کلروکین",
                           "syn": ["plaquenil"], "ind": "مالاریا و بیماری‌های خودایمنی",
                           "cls": "آمینوکینولین", "cat": "antimalarial"},

    # ---- ضدویروس بیشتر ----
    "favipiravir": {"smiles": "NC(=O)c1nc(F)cnc1O", "fa": "فاوی‌پیراویر", "syn": ["avigan"],
                    "ind": "آنفلوانزا و برخی عفونت‌های ویروسی", "cls": "مهارکننده RNA پلیمراز",
                    "cat": "antiviral"},
}

# Catalogue used for name resolution (non-drug compounds only here).
_NONDRUG = dict(_COMPOUNDS)


def _build_index() -> dict[str, str]:
    idx: dict[str, str] = {}
    for source in (_COMPOUNDS, _DRUGS):
        for canonical, info in source.items():
            smiles = info["smiles"]
            idx[canonical.lower()] = smiles
            if info.get("fa"):
                idx[info["fa"].strip().lower()] = smiles
            for syn in info.get("syn", []):
                idx[syn.strip().lower()] = smiles
    return idx


_INDEX = _build_index()


def resolve_name(text: str) -> str | None:
    """Return SMILES for a known compound/drug name, else None."""
    return _INDEX.get(text.strip().lower())


def all_compounds() -> list[dict]:
    """List the full catalogue (compounds + drugs) for the UI."""
    out = [{"name": n, "name_fa": i.get("fa"), "smiles": i["smiles"], "is_drug": False}
           for n, i in _COMPOUNDS.items()]
    out += [{"name": n, "name_fa": i.get("fa"), "smiles": i["smiles"], "is_drug": True,
             "category": i.get("cat")} for n, i in _DRUGS.items()]
    return out


def drug_info(name: str) -> dict | None:
    """Therapeutic metadata for a drug by any of its names."""
    key = name.strip().lower()
    for canonical, info in _DRUGS.items():
        names = {canonical.lower(), (info.get("fa") or "").lower()}
        names |= {s.lower() for s in info.get("syn", [])}
        if key in names:
            return {"name": canonical, "name_fa": info.get("fa"),
                    "indication_fa": info.get("ind"), "drug_class_fa": info.get("cls"),
                    "category": info.get("cat"), "smiles": info["smiles"]}
    return None


def all_drugs() -> list[dict]:
    """All drugs grouped-friendly list with therapeutic metadata."""
    return [
        {"name": n, "name_fa": i.get("fa"), "smiles": i["smiles"],
         "indication_fa": i.get("ind"), "drug_class_fa": i.get("cls"), "category": i.get("cat")}
        for n, i in _DRUGS.items()
    ]
