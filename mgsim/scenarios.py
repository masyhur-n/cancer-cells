"""Scenario construction: deterministic grid and Latin-hypercube Monte Carlo."""

import itertools

import numpy as np

from .parameters import (ALLOYS, BONE, COATINGS, CONFINEMENT, HORIZON_DAYS, IMPLANT,
                         PERFUSION_LEVELS, RESIDUAL, TISSUE, TUMOUR,
                         TUMOUR_PHENOTYPES, Range)


def _draw(rng_val, r: Range):
    """Map u in [0,1] onto a Range (log-uniform if r.log)."""
    if rng_val is None:
        return r.nominal
    if r.lo == r.hi:
        return r.lo
    if r.log:
        return float(np.exp(np.log(r.lo) + rng_val * (np.log(r.hi) - np.log(r.lo))))
    return float(r.lo + rng_val * (r.hi - r.lo))


def build_params(alloy, coating, perfusion="normal", phenotype="pH-sensitive",
                 residual="R1-marginal", form="dense", confinement="open", u=None):
    """Flat parameter dict. ``u`` maps parameter names to uniforms in [0,1]
    (Monte Carlo); missing names use the nominal value."""
    u = u or {}
    g = lambda name, r: _draw(u.get(name), r)
    co = COATINGS[coating]
    p = {
        "alloy": alloy, "coating": coating, "perfusion": perfusion,
        "phenotype": phenotype, "residual": residual, "form": form,
        "horizon": HORIZON_DAYS,
        "cr_alloy": g("cr_alloy", ALLOYS[alloy]),
        "m0": g("m0", co.m0), "m_late": g("m_late", co.m_late),
        "t_break": g("t_break", co.t_break), "width": co.width,
        "ph_cap": g("ph_cap", co.ph_cap), "bone_bonus": g("bone_bonus", co.bone_bonus),
        "half_thickness_mm": g("half_thickness_mm", IMPLANT["half_thickness_mm"]),
        "area_factor": IMPLANT["area_factor"][form],
        "T0": RESIDUAL[residual],
        "confinement": confinement,
        "confined": CONFINEMENT[confinement],
    }
    for k, r in TISSUE.items():
        p[k] = g(k, r)
    if perfusion != "sampled":
        p["k_clear"] = PERFUSION_LEVELS[perfusion]
    for k, r in TUMOUR.items():
        p[k] = g(k, r)
    for k, r in BONE.items():
        p["g0" if k == "g0_per_day" else k] = g(k, r)

    if phenotype != "sampled":
        ph = TUMOUR_PHENOTYPES[phenotype]
        p["pH50_static"] += ph["ph_shift"]
        p["pH50_kill"] += ph["ph_shift"]
        p["k_kill_max"] *= ph["kill_mult"]
        p["a_mg_pro"] = max(p["a_mg_pro"], ph["a_mg_pro"])
    return p


GRID_AXES = dict(
    alloy=list(ALLOYS),
    coating=list(COATINGS),
    perfusion=list(PERFUSION_LEVELS),
    phenotype=list(TUMOUR_PHENOTYPES),
    residual=list(RESIDUAL),
    form=list(IMPLANT["area_factor"]),
    confinement=list(CONFINEMENT),
)


def grid():
    keys = list(GRID_AXES)
    for combo in itertools.product(*GRID_AXES.values()):
        yield dict(zip(keys, combo))


# Names sampled in the Monte Carlo (each gets its own LHS column).
MC_NAMES = (["cr_alloy", "m0", "m_late", "t_break", "ph_cap", "bone_bonus",
             "half_thickness_mm"] + list(TISSUE) + list(TUMOUR) + list(BONE) + ["confined"])


def latin_hypercube(n, d, rng):
    u = (rng.permuted(np.tile(np.arange(n), (d, 1)), axis=1).T + rng.random((n, d))) / n
    return u


def monte_carlo_params(alloy, coating, n, seed=0, residual="R1-marginal", form="dense"):
    """n LHS samples for one alloy x coating. Perfusion and tumour biology are
    sampled from their full ranges (perfusion="sampled", phenotype="sampled");
    the margin is confined in half of the samples."""
    rng = np.random.default_rng(seed)
    U = latin_hypercube(n, len(MC_NAMES), rng)
    out = []
    for row in U:
        u = dict(zip(MC_NAMES, row))
        conf = "confined" if u["confined"] < 0.5 else "open"
        out.append(build_params(alloy, coating, perfusion="sampled", phenotype="sampled",
                                residual=residual, form=form, confinement=conf, u=u))
    return out
