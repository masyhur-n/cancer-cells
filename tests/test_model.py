import numpy as np

from mgsim import build_params, simulate
from mgsim.parameters import flux_to_h2_ml_cm2_day, mm_per_year_to_flux


def test_stoichiometry_units():
    # 1 mm/y of Mg ~ 2.27e-6 mol m^-2 s^-1 ~ 0.5 mL H2 cm^-2 day^-1
    j = mm_per_year_to_flux(1.0)
    assert np.isclose(j, 2.27e-6, rtol=0.01)
    assert np.isclose(flux_to_h2_ml_cm2_day(j), 0.50, rtol=0.02)


def test_inert_implant_lets_tumour_regrow():
    p = build_params("WE43", "Bare")
    p["cr_alloy"] = 0.0
    r = simulate(p)
    assert r["outcome"] == "Early recurrence"
    assert np.isclose(r["peak_interface_pH"], 7.4, atol=1e-6)


def test_more_corrosion_raises_ph_and_mg():
    lo = simulate(build_params("WE43", "PEO+HA-sealed", "low", confinement="confined"))
    hi = simulate(build_params("WE43", "Bare", "low", confinement="confined"))
    assert hi["peak_interface_pH"] > lo["peak_interface_pH"]
    assert hi["peak_interface_mg_mM"] > lo["peak_interface_mg_mM"]


def test_ha_layer_caps_interfacial_ph():
    p = build_params("ZK60", "PEO+HA", "very-low", form="porous", confinement="confined")
    p["t_break"] = 1e6        # layer never breaks
    r = simulate(p)
    assert r["peak_interface_pH"] <= p["ph_cap"] + 1e-9


def test_perfusion_dilutes_alkalinity():
    slow = simulate(build_params("AZ31", "Bare", "very-low", confinement="confined"))
    fast = simulate(build_params("AZ31", "Bare", "high", confinement="confined"))
    assert slow["peak_interface_pH"] > fast["peak_interface_pH"] + 0.5


def test_strong_alkalinity_controls_tumour_but_hurts_bone():
    p = build_params("ZK60", "Bare", "very-low", "pH-sensitive", form="porous", confinement="confined")
    p["h2_tol"] = 1e9
    p["half_thickness_mm"] = 50.0
    r = simulate(p)
    assert r["eradicated"]
    assert not r["bone_ok"]
    assert r["outcome"] == "Control, bone impaired"
