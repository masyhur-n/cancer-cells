"""Run every study and write results/ (CSV + summary.md) and figures/.

    python run_all.py            # full run (~10 min on 4 cores)
    python run_all.py --quick    # small Monte Carlo for a fast check
"""

import argparse
import csv
import os
from collections import Counter
from multiprocessing import Pool

import numpy as np

import mgsim
from mgsim.model import OUTCOMES, OUTCOME_HELP
from mgsim.parameters import ALLOYS, COATINGS, PERFUSION_LEVELS
from mgsim.scenarios import MC_NAMES, build_params, grid, monte_carlo_params
from mgsim import plots

ROOT = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")

META_KEYS = ["alloy", "coating", "perfusion", "phenotype", "residual", "form", "confinement"]
OUT_KEYS = ["outcome", "biology_outcome", "gas_failure", "structural_failure", "bone_ok", "eradicated", "recurrence_day", "log10_tumour_ratio",
            "bone_score_180", "peak_interface_pH", "mean_margin_pH", "peak_interface_mg_mM",
            "peak_h2_ml_cm2_day", "implant_remaining_84d", "implant_exhausted_day"]


def _run(p):
    r = mgsim.simulate(p)
    return {**{k: p[k] for k in META_KEYS}, **{k: r[k] for k in OUT_KEYS},
            **{f"in_{k}": p[k] for k in MC_NAMES if k in p and k != "confined"},
            "in_confined": int(p["confined"]), "in_area_factor": p["area_factor"]}


def _write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})


def _fmt_pct(counter, total, classes=OUTCOMES):
    return " | ".join(f"{100 * counter.get(o, 0) / total:.0f}%" for o in classes)


BIO = OUTCOMES[:5]


def window_sweep(pool, n=28):
    """Effective corrosion rate x perfusion for a generic, non-breaking coating."""
    crs = np.logspace(np.log10(0.003), np.log10(3.0), n)
    ks = np.logspace(np.log10(2e-6), np.log10(1e-3), n)
    jobs, idx = [], []
    for conf in ("open", "confined"):
        for i, cr in enumerate(crs):
            for j, k in enumerate(ks):
                p = build_params("WE43", "Bare", "normal", "pH-sensitive", "R1-marginal",
                                 "dense", conf)
                p.update(cr_alloy=float(cr), k_clear=float(k), half_thickness_mm=50.0,
                         h2_tol=1e9)   # isolate chemistry/biology from depletion & gas
                jobs.append(p)
                idx.append((conf, i, j))
    res = pool.map(mgsim.simulate, jobs, chunksize=8)
    grids = {c: np.empty((n, n), dtype=object) for c in ("open", "confined")}
    ph = {c: np.zeros((n, n)) for c in ("open", "confined")}
    for (c, i, j), r in zip(idx, res):
        grids[c][i, j] = r["outcome"]
        ph[c][i, j] = r["peak_interface_pH"]
    return crs, ks, grids, ph


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--mc", type=int, default=160, help="Monte Carlo samples per alloy x coating")
    args = ap.parse_args()
    n_mc = 24 if args.quick else args.mc
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)

    with Pool() as pool:
        # 1. deterministic scenario grid ---------------------------------
        jobs = [build_params(**g) for g in grid()]
        print(f"grid: {len(jobs)} scenarios")
        grid_rows = pool.map(_run, jobs, chunksize=8)
        _write(os.path.join(RES, "scenario_grid.csv"), grid_rows)

        # 2. Monte Carlo ------------------------------------------------
        mc_jobs = []
        for ai, a in enumerate(ALLOYS):
            for ci, c in enumerate(COATINGS):
                mc_jobs += monte_carlo_params(a, c, n_mc, seed=1000 + 10 * ai + ci)
                mc_jobs += monte_carlo_params(a, c, n_mc // 2, seed=5000 + 10 * ai + ci, form="porous")
        print(f"monte carlo: {len(mc_jobs)} runs")
        mc_rows = pool.map(_run, mc_jobs, chunksize=8)
        _write(os.path.join(RES, "monte_carlo.csv"), mc_rows)

        # 3. therapeutic-window map ------------------------------------
        print("window sweep")
        crs, ks, wgrids, wph = window_sweep(pool, n=14 if args.quick else 28)

    # 4. example time courses ------------------------------------------
    traces = {}
    for c in COATINGS:
        p = build_params("WE43", c, "low", "pH-sensitive", "R1-marginal", "porous", "confined")
        traces[c] = mgsim.simulate(p, record=True)["trace"]

    # figures ------------------------------------------------------------
    plots.time_courses(traces, os.path.join(FIG, "fig1_time_courses.png"))
    plots.grid_heatmap(grid_rows, os.path.join(FIG, "fig2_scenario_grid.png"))
    plots.mc_bars(mc_rows, os.path.join(FIG, "fig3_monte_carlo_outcomes.png"))
    plots.mc_bars(mc_rows, os.path.join(FIG, "fig3b_monte_carlo_tumour_bone.png"), key="biology_outcome")
    plots.window_map(crs, ks, wgrids, wph, os.path.join(FIG, "fig4_therapeutic_window.png"))
    sens = plots.sensitivity(mc_rows, os.path.join(FIG, "fig5_sensitivity.png"))

    # summary -------------------------------------------------------------
    lines = ["# Simulation results (auto-generated by `run_all.py`)", "",
             "Outcome classes:", ""]
    lines += [f"- **{o}** - {OUTCOME_HELP[o]}" for o in OUTCOMES]
    hdr = "| " + " | ".join(OUTCOMES) + " |"
    sep = "|" + "---|" * (len(OUTCOMES) + 1)

    lines += ["", f"## 1. Deterministic grid ({len(grid_rows)} scenarios)", "",
              "Share of scenarios per outcome, by coating (all alloys, perfusion levels, "
              "phenotypes, residual burdens, implant forms and confinements):", "",
              "| Coating " + hdr, sep]
    for c in COATINGS:
        sub = [r for r in grid_rows if r["coating"] == c]
        lines.append(f"| {c} | " + _fmt_pct(Counter(r["outcome"] for r in sub), len(sub)) + " |")
    lines += ["", "By perfusion level:", "", "| Perfusion " + hdr, sep]
    for pl in PERFUSION_LEVELS:
        sub = [r for r in grid_rows if r["perfusion"] == pl]
        lines.append(f"| {pl} | " + _fmt_pct(Counter(r["outcome"] for r in sub), len(sub)) + " |")
    lines += ["", "By margin confinement x implant form:", "", "| Geometry " + hdr, sep]
    for conf in ("open", "confined"):
        for form in ("dense", "porous"):
            sub = [r for r in grid_rows if r["confinement"] == conf and r["form"] == form]
            lines.append(f"| {conf}, {form} | " + _fmt_pct(Counter(r["outcome"] for r in sub), len(sub)) + " |")
    lines += ["", "By tumour phenotype:", "", "| Phenotype " + hdr, sep]
    for ph in ("pH-sensitive", "pH-resistant", "Mg-stimulated"):
        sub = [r for r in grid_rows if r["phenotype"] == ph]
        lines.append(f"| {ph} | " + _fmt_pct(Counter(r["outcome"] for r in sub), len(sub)) + " |")

    lines += ["", f"## 2. Monte Carlo ({len(mc_rows)} runs, Latin hypercube over all parameter ranges)", "",
              "Outcome probability per alloy x coating (dense + porous pooled):", "",
              "| Alloy | Coating " + hdr, "|---" + sep]
    best = []
    for a in ALLOYS:
        for c in COATINGS:
            sub = [r for r in mc_rows if r["alloy"] == a and r["coating"] == c]
            cnt = Counter(r["outcome"] for r in sub)
            best.append((cnt.get("Dual success", 0) / len(sub), a, c))
            lines.append(f"| {a} | {c} | " + _fmt_pct(cnt, len(sub)) + " |")
    best.sort(reverse=True)

    lines += ["", "### Tumour/bone outcome ignoring gas and structural limits", "",
              "The H2 criterion (Song2007-based tolerance) drives nearly all *Premature degradation* "
              "calls. The table below shows what the tumour and bone would do if the gas were "
              "tolerated or vented, plus the share of runs that also breach the gas limit:", "",
              "| Alloy | Coating | " + " | ".join(BIO) + " | gas limit breached |", "|---|---" + "|---" * (len(BIO) + 1) + "|"]
    best_bio = []
    for a in ALLOYS:
        for c in COATINGS:
            sub = [r for r in mc_rows if r["alloy"] == a and r["coating"] == c]
            cnt = Counter(r["biology_outcome"] for r in sub)
            gas = sum(r["gas_failure"] for r in sub) / len(sub)
            best_bio.append((cnt.get("Dual success", 0) / len(sub), a, c))
            lines.append(f"| {a} | {c} | " + _fmt_pct(cnt, len(sub), BIO) + f" | {100 * gas:.0f}% |")
    best_bio.sort(reverse=True)
    lines += ["", "Highest probability of tumour/bone **Dual success** (gas ignored):", ""]
    lines += [f"{i + 1}. {a} + {c}: {100 * pr:.1f}%" for i, (pr, a, c) in enumerate(best_bio[:5])]
    lines += ["", "Grid, tumour/bone outcome by perfusion x confinement (gas ignored):", "",
              "| Perfusion | Margin | " + " | ".join(BIO) + " |", "|---|---" + "|---" * len(BIO) + "|"]
    for pl in PERFUSION_LEVELS:
        for conf in ("open", "confined"):
            sub = [r for r in grid_rows if r["perfusion"] == pl and r["confinement"] == conf]
            lines.append(f"| {pl} | {conf} | " + _fmt_pct(Counter(r["biology_outcome"] for r in sub), len(sub), BIO) + " |")
    lines += ["", "Highest probability of **Dual success**:", ""]
    lines += [f"{i + 1}. {a} + {c}: {100 * pr:.1f}%" for i, (pr, a, c) in enumerate(best[:5])]

    lines += ["", "## 3. Sensitivity (Spearman rank correlation, Monte Carlo inputs)", "",
              "| Input | vs log10 tumour ratio | vs bone score (d180) |", "|---|---|---|"]
    for name, rt, rb in sens:
        lines.append(f"| {name} | {rt:+.2f} | {rb:+.2f} |")

    lines += ["", "## Figures", "",
              "- `figures/fig1_time_courses.png` - WE43 porous implant, low perfusion, confined margin, each coating",
              "- `figures/fig2_scenario_grid.png` - outcome of every grid scenario",
              "- `figures/fig3_monte_carlo_outcomes.png` - outcome probabilities per alloy x coating",
              "- `figures/fig3b_monte_carlo_tumour_bone.png` - same, tumour/bone outcome with gas limits ignored",
              "- `figures/fig4_therapeutic_window.png` - effective corrosion rate x perfusion map",
              "- `figures/fig5_sensitivity.png` - which uncertain inputs drive each outcome", ""]
    with open(os.path.join(RES, "summary.md"), "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
