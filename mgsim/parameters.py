"""Literature-informed parameter sets for the Mg-implant / residual-tumour model.

Every quantity carries either a citation key (see REFERENCES.md) or the tag
``ASSUMPTION`` when no direct measurement exists and the value is a modelling
choice. Ranges are what the Monte Carlo samples from; ``nominal`` is used for
the deterministic scenario grid.

Units: SI unless stated (concentrations in mol/m^3 == mM, time in days for the
biological rates, seconds inside the transport solver).
"""

from dataclasses import dataclass, field

# --------------------------------------------------------------------------
# Physical constants
# --------------------------------------------------------------------------
RHO_MG = 1740.0          # kg/m^3, density of Mg                      [Staiger2006]
M_MG = 0.024305          # kg/mol
V_H2_37C = 25.4e3        # mL/mol, molar gas volume at 37 C, 1 atm
SECONDS_PER_YEAR = 3.156e7
SECONDS_PER_DAY = 86400.0
PH_PHYS = 7.4
MG_SERUM = 0.8           # mM, normal serum Mg 0.7-1.05 mM          [JahnenDechent2012]


def mm_per_year_to_flux(cr_mm_y: float) -> float:
    """Corrosion penetration rate (mm/y) -> Mg dissolution flux (mol m^-2 s^-1).

    Mg -> Mg2+ + 2e-;  2H2O + 2e- -> H2 + 2OH-                       [Song1999]
    so each mol of Mg releases 1 Mg2+, 2 OH- and 1 H2.
    """
    return cr_mm_y * 1e-3 / SECONDS_PER_YEAR * RHO_MG / M_MG


def flux_to_h2_ml_cm2_day(j: float) -> float:
    return j * V_H2_37C * SECONDS_PER_DAY / 1e4


@dataclass(frozen=True)
class Range:
    lo: float
    hi: float
    nominal: float
    source: str
    log: bool = False   # sample log-uniformly


# --------------------------------------------------------------------------
# Alloys: in vivo-equivalent bare corrosion rate (mm/y)
# Wide ranges on purpose - in vivo rates are typically 1-4 orders lower than
# in vitro and depend on site and animal model             [Kirkland2012, Witte2005]
# --------------------------------------------------------------------------
ALLOYS = {
    "HP-Mg":    Range(0.2, 1.0, 0.45, "Zheng2014; Kirkland2012", log=True),
    "AZ31":     Range(0.3, 1.5, 0.70, "Witte2005; Zheng2014", log=True),
    "WE43":     Range(0.2, 1.0, 0.40, "Witte2005; Kraus2012", log=True),
    "ZK60":     Range(0.5, 2.5, 1.10, "Zheng2014", log=True),
    "Mg-Zn-Ca": Range(0.2, 1.2, 0.50, "Zheng2014; Wu2016", log=True),
}


# --------------------------------------------------------------------------
# Coatings. The coating multiplies the alloy rate by m(t):
#   m(t) = m0 + (m_late - m0) * sigmoid((t - t_break) / w)
# PEO: 3x-35x initial reduction (icorr / H2 evolution)
#      [Narayanan2014, Matykina2016, Gu2011, FattahAlhosseini2023]
# HA top-coat: further barrier + a Ca/P surface layer that caps interfacial
#      pH below ~8.5-9 while present [Lamaka2018] and is osteoconductive
#      [Hornberger2012, Song2010].
# pH_cap: ceiling of interfacial pH; Mg(OH)2 equilibrium ~10.4 for bare Mg.
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Coating:
    m0: Range
    m_late: Range
    t_break: Range          # days
    width: float            # days
    ph_cap: Range           # while intact
    bone_bonus: Range       # osteoconduction multiplier
    source: str


_CAP_MGOH2 = Range(10.2, 10.6, 10.4, "Song1999 (Mg(OH)2 equilibrium)")
_NO_BONUS = Range(1.0, 1.0, 1.0, "-")

COATINGS = {
    "Bare": Coating(
        m0=Range(1.0, 1.0, 1.0, "-"), m_late=Range(1.0, 1.0, 1.0, "-"),
        t_break=Range(0.0, 0.0, 0.0, "-"), width=1.0,
        ph_cap=_CAP_MGOH2, bone_bonus=_NO_BONUS, source="-"),
    "HA-only": Coating(
        m0=Range(0.15, 0.6, 0.35, "Song2010; Hornberger2012"),
        m_late=Range(0.8, 1.0, 0.9, "ASSUMPTION"),
        t_break=Range(5, 30, 14, "Song2010; Hornberger2012"), width=5.0,
        ph_cap=Range(8.6, 9.4, 9.0, "Lamaka2018"),
        bone_bonus=Range(1.05, 1.3, 1.15, "Hornberger2012"),
        source="electrodeposited Ca-P"),
    "PEO": Coating(
        m0=Range(0.03, 0.33, 0.10, "Matykina2016 (35x H2); Gu2011; Narayanan2014", log=True),
        m_late=Range(0.4, 0.9, 0.6, "Matykina2016 (steady state wk 6); ASSUMPTION"),
        t_break=Range(14, 60, 35, "Matykina2016; Narayanan2014"), width=7.0,
        ph_cap=_CAP_MGOH2,
        bone_bonus=Range(1.0, 1.1, 1.03, "Gu2011"),
        source="PEO / MAO"),
    "PEO+HA": Coating(
        m0=Range(0.01, 0.15, 0.04, "FattahAlhosseini2023; Narayanan2014", log=True),
        m_late=Range(0.3, 0.8, 0.5, "ASSUMPTION"),
        t_break=Range(30, 120, 60, "FattahAlhosseini2023; ASSUMPTION"), width=10.0,
        ph_cap=Range(8.5, 9.2, 8.8, "Lamaka2018"),
        bone_bonus=Range(1.1, 1.4, 1.25, "Hornberger2012; FattahAlhosseini2023"),
        source="PEO with HA incorporated / top-coat"),
    "PEO+HA-sealed": Coating(
        m0=Range(0.003, 0.05, 0.012, "Narayanan2014 (sealing strategies)", log=True),
        m_late=Range(0.2, 0.6, 0.35, "ASSUMPTION"),
        t_break=Range(90, 220, 150, "Narayanan2014; ASSUMPTION"), width=15.0,
        ph_cap=Range(8.3, 9.0, 8.6, "Lamaka2018"),
        bone_bonus=Range(1.1, 1.4, 1.25, "Hornberger2012"),
        source="PEO+HA with polymer/sol-gel sealing"),
}


# --------------------------------------------------------------------------
# Tissue / transport
# --------------------------------------------------------------------------
TISSUE = {
    # Effective diffusivities in peri-implant tissue (free water: Mg2+ 7.1e-10)
    "D_mg":   Range(2e-10, 7e-10, 4e-10, "CRC water diffusivity / tortuosity ASSUMPTION"),
    "D_base": Range(5e-10, 1.5e-9, 1e-9, "HCO3-/OH- in tissue ASSUMPTION"),
    # First-order perfusion clearance of solutes to plasma (1/s). Low end =
    # irradiated / post-resection poorly vascularised bed.
    "k_clear": Range(2e-6, 1e-3, 1e-4, "ASSUMPTION (bone perfusion order of magnitude)", log=True),
    # Bicarbonate buffer capacity, open system ~2.3*[HCO3-] ~ 55 mM/pH
    "beta":   Range(20.0, 60.0, 40.0, "Swietach2014"),
    "L":      Range(5e-3, 5e-3, 5e-3, "model margin depth"),
}

PERFUSION_LEVELS = {"very-low": 3e-6, "low": 2e-5, "normal": 1e-4, "high": 5e-4}
# Far-edge boundary condition of the margin
CONFINEMENT = {"open": False, "confined": True}


# --------------------------------------------------------------------------
# Residual tumour (normalised density T in [0,1] of carrying capacity)
# --------------------------------------------------------------------------
TUMOUR = {
    "doubling_days": Range(4.0, 30.0, 10.0, "ASSUMPTION (in vivo sarcoma regrowth)", log=True),
    "D_T":           Range(5e-15, 5e-14, 1.5e-14, "GatenbyGawlinski1996 order", log=True),
    # Cytostatic (dormancy) response to alkalinisation      [Globig2022]
    "pH50_static":   Range(7.8, 8.4, 8.05, "Globig2022"),
    "s_max":         Range(0.6, 1.0, 0.9, "Globig2022"),
    # Cytotoxic kill at strong alkalinity                    [Zhang2012, Li2014]
    "pH50_kill":     Range(8.3, 9.3, 8.7, "Zhang2012; Li2014"),
    "k_kill_max":    Range(0.3, 3.0, 1.0, "ASSUMPTION calibrated to Zhang2012/Li2014", log=True),
    # Mg2+ cytotoxicity (> ~10 mM reduces osteosarcoma cell number) [Burmester2014]
    "IC50_mg":       Range(10.0, 40.0, 20.0, "Burmester2014; Zhang2012"),
    "k_mg_max":      Range(0.05, 0.6, 0.2, "ASSUMPTION", log=True),
    # Mg2+ as a proliferation co-factor                      [WolfTrapani2008, CastiglioniMaier2011]
    "a_mg_pro":      Range(0.0, 0.3, 0.0, "WolfTrapani2008; CastiglioniMaier2011"),
    # Tumour acidosis: pHe depression at full density       [Webb2011, GatenbyGawlinski1996]
    "dpH_acid":      Range(0.3, 0.7, 0.5, "Webb2011; Swietach2014"),
}

TUMOUR_PHENOTYPES = {
    # (multiplier on pH50 offsets, k_kill_max multiplier, a_mg_pro)
    "pH-sensitive": dict(ph_shift=-0.15, kill_mult=1.5, a_mg_pro=0.0),
    "pH-resistant": dict(ph_shift=+0.35, kill_mult=0.4, a_mg_pro=0.0),
    "Mg-stimulated": dict(ph_shift=+0.10, kill_mult=0.8, a_mg_pro=0.3),
}

RESIDUAL = {
    # Mean initial residual density in the first 1 mm of margin
    "R0-micro":   1e-5,
    "R1-marginal": 1e-3,
}


# --------------------------------------------------------------------------
# Bone formation (osteogenic index, relative to an inert Ti implant)
# --------------------------------------------------------------------------
BONE = {
    "g0_per_day":   Range(0.008, 0.02, 0.012, "ASSUMPTION (Ti reaches ~88% by d180)"),
    # Mg2+ stimulation 2-10 mM                [Yoshizawa2014, Zhang2016]
    "mg_stim_amp":  Range(0.2, 0.8, 0.5, "Yoshizawa2014; Zhang2016"),
    "mg_stim_K":    Range(1.0, 5.0, 2.5, "Yoshizawa2014"),
    # Mg2+ inhibition at high concentration    [Wang2015, Burmester2014]
    "mg_tox_IC50":  Range(15.0, 40.0, 25.0, "Wang2015; Burmester2014"),
    # Alkaline optimum ~pH 8.0-8.4 in vitro    [Galow2017]
    "ph_opt":       Range(7.9, 8.4, 8.2, "Galow2017"),
    "ph_amp":       Range(0.1, 0.6, 0.3, "Galow2017 (in vivo damped, ASSUMPTION)"),
    "ph_tox50":     Range(8.9, 9.8, 9.3, "Galow2017 (>9 untested); ASSUMPTION"),
    # H2 tolerance, 0.01 mL/cm2/day widely cited lower bound [Song2007]
    # Upper end reflects clinical Mg screws that show gas radiolucency yet heal
    # [Witte2010]; gas failure = >5x tolerance sustained >= 28 days.
    "h2_tol":       Range(0.01, 0.3, 0.06, "Song2007 (lower); Witte2010; ASSUMPTION", log=True),
}

IMPLANT = {
    "half_thickness_mm": Range(1.0, 2.5, 1.5, "ASSUMPTION (plate/scaffold strut)"),
    "area_factor":       {"dense": 1.0, "porous": 3.0},   # Hakimi2016 porous EW62
}

HORIZON_DAYS = 365.0
