"""1-D reaction-diffusion model of a degrading Mg implant in a resection margin.

Domain: x in [0, L] measured from the implant surface into the peri-implant
tissue that may contain residual tumour cells after resection. The far edge is
either "open" (healthy, perfused tissue holds serum levels) or "confined"
(zero flux: fibrous capsule, closed medullary cavity, haematoma).

State per grid cell
  c  excess Mg2+ over serum level (mM)
  b  excess base (mM of OH- equivalents; negative = acid). pH = 7.4 + b/beta
  T  residual tumour density, normalised to carrying capacity (0..1)
  O  osteogenic / new-bone index (0..1)

Chemistry: each mol of corroded Mg releases 1 Mg2+, 2 OH-, 1 H2 [Song1999].
Solutes diffuse, are cleared by perfusion (first order to plasma) and the
interfacial pH is capped by precipitation of Mg(OH)2 or, while an HA layer
persists, a Ca-P layer [Lamaka2018]. Tumour cells produce acid
[GatenbyGawlinski1996, Webb2011], which counteracts alkalinisation.

Biology:
  tumour: logistic growth, cytostatic (dormancy) response at moderate
  alkalinity [Globig2022], cytotoxic kill at strong alkalinity [Zhang2012,
  Li2014], Mg2+ cytotoxicity above ~10 mM [Burmester2014], optional Mg2+
  pro-proliferative effect [WolfTrapani2008].
  bone: osteogenic rate modulated by Mg2+ (biphasic) [Yoshizawa2014,
  Zhang2016, Wang2015], pH (optimum ~8.2) [Galow2017], H2 gas near the
  implant [Song2007], and the HA osteoconduction bonus [Hornberger2012];
  space occupied by tumour is unavailable for bone.

Numerics: backward-Euler diffusion (pre-inverted matrices), operator-split
exponential integration of the reaction terms, dt = 0.2 day.
"""

import numpy as np

from .parameters import (MG_SERUM, PH_PHYS, SECONDS_PER_DAY,
                         flux_to_h2_ml_cm2_day, mm_per_year_to_flux)

N_CELLS = 50
DT_DAYS = 0.2
MGOH2_CAP = 10.4

OUTCOMES = [
    "Dual success",
    "Dormant control",
    "Control, bone impaired",
    "Late recurrence",
    "Early recurrence",
    "Premature degradation",
]

OUTCOME_HELP = {
    "Dual success": "residual tumour eradicated and bone formation >= 80% of a Ti control",
    "Dormant control": "tumour suppressed but viable cells remain (relapse risk once Mg is spent); bone OK",
    "Control, bone impaired": "tumour eradicated or suppressed, but alkalinity / Mg2+ / H2 cut bone formation < 80%",
    "Late recurrence": "detectable regrowth between day 180 and 365 (after coating breakdown or implant exhaustion)",
    "Early recurrence": "detectable regrowth within 180 days (degradation too slow / tumour too resistant)",
    "Premature degradation": "<50% implant left at 12 weeks or H2 > 5x tolerance for >= 28 d (structural / gas failure)",
}


def _hill(x, k, n):
    x = np.maximum(x, 0.0)
    return x ** n / (x ** n + k ** n)


def _sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -50, 50)))


def _implicit_inverse(n, dx, D, k, dt, right_dirichlet):
    """Inverse of (I + dt*(k - D*Lap)) for a finite-volume grid.

    Left face: zero diffusive flux (sources added explicitly).
    Right face: Dirichlet zero excess (healthy perfused tissue) or zero flux.
    """
    lap = np.zeros((n, n))
    for i in range(n):
        if i > 0:
            lap[i, i - 1] = 1.0
            lap[i, i] -= 1.0
        if i < n - 1:
            lap[i, i + 1] = 1.0
            lap[i, i] -= 1.0
    if right_dirichlet:
        lap[n - 1, n - 1] -= 2.0     # ghost node = -c_{N-1} -> face value 0
    lap /= dx * dx
    a = np.eye(n) * (1.0 + dt * k) - dt * D * lap
    return np.linalg.inv(a)


def simulate(p, record=False):
    """Run one simulation. ``p`` is a flat dict (see scenarios.build_params)."""
    L = p["L"]
    n = N_CELLS
    dx = L / n
    x = (np.arange(n) + 0.5) * dx
    dt_d = DT_DAYS
    dt_s = dt_d * SECONDS_PER_DAY
    steps = int(round(p["horizon"] / dt_d))

    open_edge = not p["confined"]
    inv_mg = _implicit_inverse(n, dx, p["D_mg"], p["k_clear"], dt_s, open_edge)
    inv_b = _implicit_inverse(n, dx, p["D_base"], p["k_clear"], dt_s, open_edge)
    inv_t = _implicit_inverse(n, dx, p["D_T"], 0.0, dt_s, False)

    beta = p["beta"]
    # Tumour acid production chosen so a fully packed, well-mixed tumour sits
    # dpH_acid below 7.4 at steady state (q / k_clear = dpH * beta).
    q_acid = p["dpH_acid"] * beta * p["k_clear"]          # mM/s at T = 1

    c = np.zeros(n)
    b = np.zeros(n)
    prof = np.exp(-x / 0.5e-3)
    prof /= prof[x <= 1e-3].mean()
    T = p["T0"] * prof
    O = np.zeros(n)
    near = np.exp(-x / 0.5e-3)            # H2 gas acts near the interface

    r = np.log(2.0) / p["doubling_days"]
    pH50s = p["pH50_static"]
    pH50k = p["pH50_kill"]
    N0 = T.sum()
    m2 = x <= 2e-3

    depth = 0.0
    j0 = mm_per_year_to_flux(p["cr_alloy"]) * p["area_factor"]

    peak_pH = PH_PHYS
    peak_mg = MG_SERUM
    peak_h2 = 0.0
    gas_over_days = 0.0
    recur_day = None
    remaining_84 = None
    exhausted_day = None
    pH_integral = 0.0
    T_max_seen = T.max()
    bone_180 = None
    rec = {k: [] for k in ("t", "pH0", "mg0", "tumour", "bone", "h2", "remaining")} if record else None

    f_ph7 = 1.0 + p["ph_amp"] * np.exp(-((PH_PHYS - p["ph_opt"]) / 0.4) ** 2)

    for s in range(1, steps + 1):
        t = s * dt_d
        sig = _sigmoid((t - p["t_break"]) / p["width"]) if p["t_break"] > 0 else 1.0
        m = p["m0"] + (p["m_late"] - p["m0"]) * sig
        cap = p["ph_cap"] + (MGOH2_CAP - p["ph_cap"]) * sig

        if depth < p["half_thickness_mm"]:
            cr_eff = p["cr_alloy"] * m                # mm/y
            depth += cr_eff * dt_d / 365.25
            j = j0 * m
        else:
            j = 0.0
            if exhausted_day is None:
                exhausted_day = t
        h2 = flux_to_h2_ml_cm2_day(j)
        peak_h2 = max(peak_h2, h2)
        if h2 > 5.0 * p["h2_tol"]:
            gas_over_days += dt_d
        if remaining_84 is None and t >= 84.0:
            remaining_84 = max(0.0, 1.0 - depth / p["half_thickness_mm"])

        # --- transport (mM; 1 mol/m^3 == 1 mM) ---
        src_mg = np.zeros(n)
        src_mg[0] = j / dx * dt_s
        src_b = -q_acid * T * dt_s
        src_b[0] += 2.0 * j / dx * dt_s
        c = inv_mg @ (c + src_mg)
        b = inv_b @ (b + src_b)

        # precipitation at the pH ceiling (Mg(OH)2 or Ca-P layer)
        b_cap = beta * (cap - PH_PHYS)
        over = np.maximum(b - b_cap, 0.0)
        b -= over
        c = np.maximum(c - 0.5 * over, 0.0)

        pH = np.clip(PH_PHYS + b / beta, 6.0, cap)
        C = MG_SERUM + c
        peak_pH = max(peak_pH, pH[0])
        peak_mg = max(peak_mg, C[0])
        pH_integral += pH[m2].mean() * dt_d

        # --- tumour ---
        T = inv_t @ T
        static = p["s_max"] * _hill(pH - PH_PHYS, pH50s - PH_PHYS, 4)
        growth = r * (1.0 + p["a_mg_pro"] * _hill(c, 3.0, 2)) * (1.0 - static)
        kill = (p["k_kill_max"] * _hill(pH - PH_PHYS, pH50k - PH_PHYS, 6)
                + p["k_mg_max"] * _hill(C, p["IC50_mg"], 2))
        T = T * np.exp((growth * (1.0 - T) - kill) * dt_d)
        T = np.clip(T, 0.0, 1.0)
        T_max_seen = max(T_max_seen, T.max())
        if recur_day is None and T[m2].mean() > 0.05:
            recur_day = t

        # --- bone ---
        f_mg = (1.0 + p["mg_stim_amp"] * _hill(c, p["mg_stim_K"], 2)) * (1.0 - _hill(C, p["mg_tox_IC50"], 3))
        f_ph = (1.0 + p["ph_amp"] * np.exp(-((pH - p["ph_opt"]) / 0.4) ** 2)) / f_ph7
        f_ph *= 1.0 - _hill(pH - PH_PHYS, p["ph_tox50"] - PH_PHYS, 8)
        f_ph *= np.where(pH < PH_PHYS, np.exp(-2.0 * (PH_PHYS - pH)), 1.0)   # acidosis
        f_gas = 1.0 - near * (1.0 - 1.0 / (1.0 + (h2 / p["h2_tol"]) ** 2))
        rate = p["g0"] * p["bone_bonus"] * f_mg * f_ph * f_gas * (1.0 - T)
        O = 1.0 - (1.0 - O) * np.exp(-np.maximum(rate, 0.0) * dt_d)
        # alkaline / hyper-Mg injury of newly formed tissue (ASSUMPTION)
        injury = 0.03 * np.maximum(_hill(pH - PH_PHYS, p["ph_tox50"] - PH_PHYS + 0.3, 8),
                                   _hill(C, 2.0 * p["mg_tox_IC50"], 3))
        O = O * np.exp(-injury * dt_d)

        if bone_180 is None and t >= 180.0:
            bone_180 = O[m2].mean() / (1.0 - np.exp(-p["g0"] * t))

        if record and (s % 5 == 0):
            rec["t"].append(t)
            rec["pH0"].append(pH[0])
            rec["mg0"].append(C[0])
            rec["tumour"].append(T.sum() / N0)
            rec["bone"].append(O[m2].mean() / (1.0 - np.exp(-p["g0"] * t)))
            rec["h2"].append(h2)
            rec["remaining"].append(max(0.0, 1.0 - depth / p["half_thickness_mm"]))

    tumour_ratio = T.sum() / N0
    eradicated = recur_day is None and T.max() < 1e-9
    bone_ok = bone_180 >= 0.8
    premature = (remaining_84 is not None and remaining_84 < 0.5) or gas_over_days >= 28.0

    if premature:
        outcome = "Premature degradation"
    elif recur_day is not None and recur_day <= 180.0:
        outcome = "Early recurrence"
    elif recur_day is not None:
        outcome = "Late recurrence"
    elif not bone_ok:
        outcome = "Control, bone impaired"
    elif eradicated:
        outcome = "Dual success"
    else:
        outcome = "Dormant control"

    out = dict(
        outcome=outcome,
        bone_ok=bool(bone_ok),
        eradicated=bool(eradicated),
        recurrence_day=recur_day if recur_day is not None else np.nan,
        log10_tumour_ratio=float(np.log10(max(tumour_ratio, 1e-300))),
        bone_score_180=float(bone_180),
        peak_interface_pH=float(peak_pH),
        mean_margin_pH=float(pH_integral / p["horizon"]),
        peak_interface_mg_mM=float(peak_mg),
        peak_h2_ml_cm2_day=float(peak_h2),
        implant_remaining_84d=float(remaining_84 if remaining_84 is not None else np.nan),
        implant_exhausted_day=exhausted_day if exhausted_day is not None else np.nan,
    )
    if record:
        out["trace"] = {k: np.asarray(v) for k, v in rec.items()}
    return out
