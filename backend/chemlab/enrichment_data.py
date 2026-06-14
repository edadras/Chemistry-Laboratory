"""
enrichment_data.py — Curated element data that augments RDKit.

`ELEMENT_ENRICHMENT[Z]` holds a dict of extra properties for element with
atomic number Z. Names (English + Persian) are provided for all 118 elements;
physical/chemical detail is curated for the elements most relevant to a
general and pharmaceutical chemistry lab.
"""
from __future__ import annotations

# Persian labels for element categories.
CATEGORY_FA = {
    "nonmetal": "نافلز",
    "noble_gas": "گاز نجیب",
    "alkali_metal": "فلز قلیایی",
    "alkaline_earth_metal": "فلز قلیایی خاکی",
    "metalloid": "شبه‌فلز",
    "halogen": "هالوژن",
    "post_transition_metal": "فلز پسواسطه",
    "transition_metal": "فلز واسطه",
    "lanthanide": "لانتانید",
    "actinide": "اکتینید",
    "unknown": "نامشخص",
}

# (English name, Persian name) for every element 1..118.
_NAMES: dict[int, tuple[str, str]] = {
    1: ("Hydrogen", "هیدروژن"), 2: ("Helium", "هلیوم"), 3: ("Lithium", "لیتیوم"),
    4: ("Beryllium", "بریلیوم"), 5: ("Boron", "بور"), 6: ("Carbon", "کربن"),
    7: ("Nitrogen", "نیتروژن"), 8: ("Oxygen", "اکسیژن"), 9: ("Fluorine", "فلوئور"),
    10: ("Neon", "نئون"), 11: ("Sodium", "سدیم"), 12: ("Magnesium", "منیزیم"),
    13: ("Aluminium", "آلومینیوم"), 14: ("Silicon", "سیلیسیم"), 15: ("Phosphorus", "فسفر"),
    16: ("Sulfur", "گوگرد"), 17: ("Chlorine", "کلر"), 18: ("Argon", "آرگون"),
    19: ("Potassium", "پتاسیم"), 20: ("Calcium", "کلسیم"), 21: ("Scandium", "اسکاندیم"),
    22: ("Titanium", "تیتانیوم"), 23: ("Vanadium", "وانادیوم"), 24: ("Chromium", "کروم"),
    25: ("Manganese", "منگنز"), 26: ("Iron", "آهن"), 27: ("Cobalt", "کبالت"),
    28: ("Nickel", "نیکل"), 29: ("Copper", "مس"), 30: ("Zinc", "روی"),
    31: ("Gallium", "گالیم"), 32: ("Germanium", "ژرمانیوم"), 33: ("Arsenic", "آرسنیک"),
    34: ("Selenium", "سلنیوم"), 35: ("Bromine", "برم"), 36: ("Krypton", "کریپتون"),
    37: ("Rubidium", "روبیدیم"), 38: ("Strontium", "استرانسیم"), 39: ("Yttrium", "ایتریم"),
    40: ("Zirconium", "زیرکونیم"), 41: ("Niobium", "نیوبیم"), 42: ("Molybdenum", "مولیبدن"),
    43: ("Technetium", "تکنسیم"), 44: ("Ruthenium", "روتنیم"), 45: ("Rhodium", "رودیم"),
    46: ("Palladium", "پالادیم"), 47: ("Silver", "نقره"), 48: ("Cadmium", "کادمیم"),
    49: ("Indium", "ایندیم"), 50: ("Tin", "قلع"), 51: ("Antimony", "آنتیموان"),
    52: ("Tellurium", "تلوریم"), 53: ("Iodine", "ید"), 54: ("Xenon", "زنون"),
    55: ("Caesium", "سزیم"), 56: ("Barium", "باریم"), 57: ("Lanthanum", "لانتانیم"),
    58: ("Cerium", "سریم"), 59: ("Praseodymium", "پرازئودیمیم"), 60: ("Neodymium", "نئودیمیم"),
    61: ("Promethium", "پرومتیم"), 62: ("Samarium", "ساماریم"), 63: ("Europium", "اوروپیم"),
    64: ("Gadolinium", "گادولینیم"), 65: ("Terbium", "تربیم"), 66: ("Dysprosium", "دیسپروزیم"),
    67: ("Holmium", "هولمیم"), 68: ("Erbium", "اربیم"), 69: ("Thulium", "تولیم"),
    70: ("Ytterbium", "ایتربیم"), 71: ("Lutetium", "لوتسیم"), 72: ("Hafnium", "هافنیم"),
    73: ("Tantalum", "تانتالیم"), 74: ("Tungsten", "تنگستن"), 75: ("Rhenium", "رنیم"),
    76: ("Osmium", "اسمیم"), 77: ("Iridium", "ایریدیم"), 78: ("Platinum", "پلاتین"),
    79: ("Gold", "طلا"), 80: ("Mercury", "جیوه"), 81: ("Thallium", "تالیم"),
    82: ("Lead", "سرب"), 83: ("Bismuth", "بیسموت"), 84: ("Polonium", "پولونیم"),
    85: ("Astatine", "استاتین"), 86: ("Radon", "رادون"), 87: ("Francium", "فرانسیم"),
    88: ("Radium", "رادیم"), 89: ("Actinium", "اکتینیم"), 90: ("Thorium", "توریم"),
    91: ("Protactinium", "پروتاکتینیم"), 92: ("Uranium", "اورانیم"), 93: ("Neptunium", "نپتونیم"),
    94: ("Plutonium", "پلوتونیم"), 95: ("Americium", "آمریسیم"), 96: ("Curium", "کوریم"),
    97: ("Berkelium", "برکلیم"), 98: ("Californium", "کالیفرنیم"), 99: ("Einsteinium", "اینشتینیم"),
    100: ("Fermium", "فرمیم"), 101: ("Mendelevium", "مندلیویم"), 102: ("Nobelium", "نوبلیم"),
    103: ("Lawrencium", "لارنسیم"), 104: ("Rutherfordium", "رادرفوردیم"), 105: ("Dubnium", "دوبنیم"),
    106: ("Seaborgium", "سیبورگیم"), 107: ("Bohrium", "بوریم"), 108: ("Hassium", "هاسیم"),
    109: ("Meitnerium", "مایتنریم"), 110: ("Darmstadtium", "دارمشتاتیم"), 111: ("Roentgenium", "رونتگنیم"),
    112: ("Copernicium", "کوپرنیسیم"), 113: ("Nihonium", "نیهونیم"), 114: ("Flerovium", "فلروویم"),
    115: ("Moscovium", "مسکوویم"), 116: ("Livermorium", "لیورموریم"), 117: ("Tennessine", "تنسین"),
    118: ("Oganesson", "اوگانسون"),
}

# Rich curated data for the elements most relevant to a chemistry/pharma lab.
# Keys: category, group, period, electron_config, electronegativity (Pauling),
# melting_point_k, boiling_point_k, density, oxidation_states, phase_stp,
# discovered, biological_role.
_RICH: dict[int, dict] = {
    1: dict(category="nonmetal", group=1, period=1, electron_config="1s1",
            electronegativity=2.20, melting_point_k=14.01, boiling_point_k=20.28,
            density=0.00008988, oxidation_states=[-1, 1], phase_stp="gas",
            discovered="1766", biological_role="جزء آب و تمام مولکول‌های آلی"),
    2: dict(category="noble_gas", group=18, period=1, electron_config="1s2",
            electronegativity=None, melting_point_k=0.95, boiling_point_k=4.22,
            density=0.0001785, oxidation_states=[0], phase_stp="gas",
            discovered="1868", biological_role="بی‌اثر"),
    6: dict(category="nonmetal", group=14, period=2, electron_config="1s2 2s2 2p2",
            electronegativity=2.55, melting_point_k=3823, boiling_point_k=4098,
            density=2.267, oxidation_states=[-4, -3, -2, -1, 1, 2, 3, 4], phase_stp="solid",
            discovered="باستان", biological_role="ستون فقرات شیمی آلی و حیات"),
    7: dict(category="nonmetal", group=15, period=2, electron_config="1s2 2s2 2p3",
            electronegativity=3.04, melting_point_k=63.15, boiling_point_k=77.36,
            density=0.0012506, oxidation_states=[-3, -2, -1, 1, 2, 3, 4, 5], phase_stp="gas",
            discovered="1772", biological_role="جزء پروتئین‌ها، DNA و آمینواسیدها"),
    8: dict(category="nonmetal", group=16, period=2, electron_config="1s2 2s2 2p4",
            electronegativity=3.44, melting_point_k=54.36, boiling_point_k=90.20,
            density=0.001429, oxidation_states=[-2, -1, 1, 2], phase_stp="gas",
            discovered="1771", biological_role="تنفس سلولی، جزء آب"),
    9: dict(category="halogen", group=17, period=2, electron_config="1s2 2s2 2p5",
            electronegativity=3.98, melting_point_k=53.53, boiling_point_k=85.03,
            density=0.001696, oxidation_states=[-1], phase_stp="gas",
            discovered="1886", biological_role="استحکام مینای دندان"),
    11: dict(category="alkali_metal", group=1, period=3, electron_config="[Ne] 3s1",
             electronegativity=0.93, melting_point_k=370.87, boiling_point_k=1156,
             density=0.971, oxidation_states=[1], phase_stp="solid",
             discovered="1807", biological_role="تعادل الکترولیت و انتقال عصبی"),
    12: dict(category="alkaline_earth_metal", group=2, period=3, electron_config="[Ne] 3s2",
             electronegativity=1.31, melting_point_k=923, boiling_point_k=1363,
             density=1.738, oxidation_states=[2], phase_stp="solid",
             discovered="1755", biological_role="کوفاکتور آنزیمی، جزء کلروفیل"),
    15: dict(category="nonmetal", group=15, period=3, electron_config="[Ne] 3s2 3p3",
             electronegativity=2.19, melting_point_k=317.30, boiling_point_k=550,
             density=1.823, oxidation_states=[-3, 3, 5], phase_stp="solid",
             discovered="1669", biological_role="جزء DNA, ATP و استخوان"),
    16: dict(category="nonmetal", group=16, period=3, electron_config="[Ne] 3s2 3p4",
             electronegativity=2.58, melting_point_k=388.36, boiling_point_k=717.87,
             density=2.067, oxidation_states=[-2, 2, 4, 6], phase_stp="solid",
             discovered="باستان", biological_role="جزء آمینواسیدهای سیستئین و متیونین"),
    17: dict(category="halogen", group=17, period=3, electron_config="[Ne] 3s2 3p5",
             electronegativity=3.16, melting_point_k=171.6, boiling_point_k=239.11,
             density=0.003214, oxidation_states=[-1, 1, 3, 5, 7], phase_stp="gas",
             discovered="1774", biological_role="تعادل الکترولیت، اسید معده (HCl)"),
    19: dict(category="alkali_metal", group=1, period=4, electron_config="[Ar] 4s1",
             electronegativity=0.82, melting_point_k=336.53, boiling_point_k=1032,
             density=0.862, oxidation_states=[1], phase_stp="solid",
             discovered="1807", biological_role="انتقال عصبی و عملکرد قلب"),
    20: dict(category="alkaline_earth_metal", group=2, period=4, electron_config="[Ar] 4s2",
             electronegativity=1.00, melting_point_k=1115, boiling_point_k=1757,
             density=1.54, oxidation_states=[2], phase_stp="solid",
             discovered="1808", biological_role="استخوان، دندان و انعقاد خون"),
    26: dict(category="transition_metal", group=8, period=4, electron_config="[Ar] 3d6 4s2",
             electronegativity=1.83, melting_point_k=1811, boiling_point_k=3134,
             density=7.874, oxidation_states=[2, 3, 6], phase_stp="solid",
             discovered="باستان", biological_role="حمل اکسیژن در هموگلوبین"),
    29: dict(category="transition_metal", group=11, period=4, electron_config="[Ar] 3d10 4s1",
             electronegativity=1.90, melting_point_k=1357.77, boiling_point_k=2835,
             density=8.96, oxidation_states=[1, 2], phase_stp="solid",
             discovered="باستان", biological_role="کوفاکتور آنزیم‌های اکسیداز"),
    30: dict(category="transition_metal", group=12, period=4, electron_config="[Ar] 3d10 4s2",
             electronegativity=1.65, melting_point_k=692.68, boiling_point_k=1180,
             density=7.134, oxidation_states=[2], phase_stp="solid",
             discovered="باستان", biological_role="کوفاکتور بیش از ۳۰۰ آنزیم"),
    35: dict(category="halogen", group=17, period=4, electron_config="[Ar] 3d10 4s2 4p5",
             electronegativity=2.96, melting_point_k=265.8, boiling_point_k=332.0,
             density=3.122, oxidation_states=[-1, 1, 3, 5, 7], phase_stp="liquid",
             discovered="1826", biological_role="کمیاب در بدن"),
    53: dict(category="halogen", group=17, period=5, electron_config="[Kr] 4d10 5s2 5p5",
             electronegativity=2.66, melting_point_k=386.85, boiling_point_k=457.4,
             density=4.93, oxidation_states=[-1, 1, 3, 5, 7], phase_stp="solid",
             discovered="1811", biological_role="ساخت هورمون‌های تیروئید"),
    79: dict(category="transition_metal", group=11, period=6, electron_config="[Xe] 4f14 5d10 6s1",
             electronegativity=2.54, melting_point_k=1337.33, boiling_point_k=3129,
             density=19.282, oxidation_states=[1, 3], phase_stp="solid",
             discovered="باستان", biological_role="بی‌اثر؛ کاربرد دارویی محدود"),
}


def _merge() -> dict[int, dict]:
    out: dict[int, dict] = {}
    for z, (name, name_fa) in _NAMES.items():
        d = {"name": name, "name_fa": name_fa}
        d.update(_RICH.get(z, {}))
        out[z] = d
    return out


ELEMENT_ENRICHMENT: dict[int, dict] = _merge()
