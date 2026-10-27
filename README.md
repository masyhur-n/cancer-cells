# Mg-alloy / PEO / HA implants after bone-tumour resection: an outcome simulator

After a bone tumour is resected, the implant that fills the defect has two
jobs that can pull against each other:

1. **Support new bone** (osteoconduction, mechanical support while the bone heals).
2. **Stop residual tumour cells at the margin from regrowing** (local recurrence).

Biodegradable Mg alloys are candidates for both jobs. Corroding Mg releases
**Mg²⁺** (osteogenic at a few mM [Yoshizawa2014, Zhang2016]) and **OH⁻**
(local alkalinisation that is cytotoxic or cytostatic to osteosarcoma cells
in vitro [Zhang2012, Li2014, Globig2022]). PEO/MAO and HA coatings slow that
release and change the interface chemistry [Narayanan2014, Lamaka2018].
Nobody has quantified how much alkalinity and Mg²⁺ actually reach residual
tumour cells in vivo. This repository is a mechanistic, literature-parameterised
model that **enumerates the possible outcomes** and shows which conditions lead
to each.

> **Scope.** This is a hypothesis-generating model. Many parameters come
> from in vitro studies or are explicit assumptions (tagged `ASSUMPTION` in
> `mgsim/parameters.py`). Read its outputs as a map of plausible regimes and of
> the measurements that would decide between them, not as predictions for a
> patient. See [Limitations](#limitations).

All citations: [REFERENCES.md](REFERENCES.md).

---

## Model

1-D reaction–diffusion along a 5 mm resection margin, measured outward from
the implant surface (`mgsim/model.py`).

```
 implant | interface |------------- peri-implant margin (0..5 mm) -------------| edge
  Mg --> Mg2+ + 2OH- + H2      diffusion + perfusion clearance             open: serum levels
  coating factor m(t)          tumour cells T(x,t)  (acid producers)       confined: zero flux
  pH ceiling (Mg(OH)2 / CaP)   new bone O(x,t)
```

| Layer | What is modelled | Basis |
|---|---|---|
| Corrosion | Flux `J = CR·ρ/M · area · m(t)`; 1 Mg²⁺ + 2 OH⁻ + 1 H₂ per Mg | [Song1999] |
| Coatings | `m(t)` goes sigmoidally from a barrier value `m0` to `m_late` at breakdown time `t_break` | PEO [Matykina2016, Narayanan2014, Gu2011]; HA [Hornberger2012, Song2010, FattahAlhosseini2023] |
| Interface pH ceiling | 10.4 (Mg(OH)₂); 8.3–9.4 while a Ca-P/HA layer persists | [Lamaka2018] |
| Transport | Implicit diffusion of Mg²⁺ and base; first-order perfusion clearance; bicarbonate buffer `pH = 7.4 + b/β` | [Swietach2014] |
| Tumour | Logistic growth; cytostatic (dormancy) effect of alkalinity; cytotoxic kill at high pH; Mg²⁺ kill above ~10 mM; optional Mg²⁺-driven proliferation; tumour cells acidify their surroundings | [Globig2022, Zhang2012, Li2014, Burmester2014, WolfTrapani2008, Webb2011, GatenbyGawlinski1996] |
| Bone | Osteogenic rate × Mg²⁺ factor (rises then falls with dose) × pH factor (optimum ~8.2, toxic above ~9.3) × H₂-gas factor × HA bonus × space not taken by tumour; alkaline injury term | [Yoshizawa2014, Zhang2016, Wang2015, Galow2017, Song2007, Hornberger2012] |
| Implant | Penetration depth vs strut half-thickness (structural reserve) | — |

### Outcome classes

Each simulation (365 days) ends in exactly one class:

| | Outcome | Definition |
|---|---|---|
| A | **Dual success** | Residual tumour eradicated, and bone formation at day 180 is at least 80% of an inert Ti control |
| B | **Dormant control** | Tumour suppressed but viable cells remain, so it can relapse once the Mg is used up; bone OK |
| C | **Control, bone impaired** | Tumour controlled, but alkalinity, Mg²⁺ or H₂ pushed bone formation below 80% |
| D | **Late recurrence** | Detectable regrowth between day 180 and 365, after coating breakdown or implant exhaustion |
| E | **Early recurrence** | Detectable regrowth within 180 days: degradation too slow, or tumour too resistant |
| F | **Premature degradation** | Less than 50% of the implant left at 12 weeks, or H₂ above 5× tolerance for at least 28 days |

### Scenario axes

| Axis | Levels |
|---|---|
| Alloy (in vivo bare rate, mm/y) | HP-Mg, AZ31, WE43, ZK60, Mg-Zn-Ca |
| Coating | Bare, HA-only, PEO, PEO+HA, PEO+HA-sealed |
| Perfusion clearance | very-low (irradiated / poorly vascular bed), low, normal, high |
| Margin geometry | open (healthy perfused tissue at 5 mm) vs confined (capsule, closed cavity, haematoma) |
| Implant form | dense vs porous scaffold (3× surface area) [Hakimi2016] |
| Tumour phenotype | pH-sensitive, pH-resistant, Mg-stimulated |
| Residual burden | R0-micro (1e-5 of carrying capacity), R1-marginal (1e-3) |

The **deterministic grid** covers all 2,400 combinations at nominal values.
The **Monte Carlo** draws a Latin hypercube over every uncertain parameter
range (25 alloy × coating pairs × 240 samples = 6,000 runs).

---

## Run it

```bash
pip install -r requirements.txt
python run_all.py          # full study, ~15-20 min on 4 cores
python run_all.py --quick  # smaller Monte Carlo for a fast check
python -m pytest -q        # sanity tests (stoichiometry, buffering, pH caps, trade-off)
```

Outputs:

- `results/scenario_grid.csv`: one row per grid scenario, with outcome, peak pH and Mg²⁺, H₂, bone score, recurrence day
- `results/monte_carlo.csv`: every Monte Carlo sample, with its sampled inputs (`in_*`)
- `results/summary.md`: auto-generated tables (outcome shares, probabilities, sensitivity)
- `figures/fig1..fig5*.png`

Single run from Python:

```python
from mgsim import build_params, simulate
p = build_params("WE43", "PEO+HA", perfusion="low", phenotype="pH-sensitive",
                 residual="R1-marginal", form="porous", confinement="confined")
r = simulate(p, record=True)
print(r["outcome"], r["peak_interface_pH"], r["bone_score_180"])
```

---

## Results

<!-- RESULTS -->

---

## Limitations

- **1-D, continuum, deterministic biology.** No immune response, angiogenesis,
  osteoclasts, mechanics or individual-cell stochastic extinction. "Eradication"
  means the continuum density fell below 1e-9 of carrying capacity.
- **Most dose–response curves are in vitro** (cell lines in HEPES- or
  extract-based media). In vivo buffering, protein adsorption and the immune
  system are likely to shift them.
- **In vivo corrosion rates and perfusion of a post-resection bed are poorly
  constrained.** That is exactly why the model sweeps them, and why they
  dominate the sensitivity analysis.
- **Coating behaviour is simplified** to a barrier factor with a single
  breakdown time. Galvanic effects, coating cracks, localised or pitting
  corrosion and alloying-element ions (Zn, rare earths, Al) are not modelled.
- **No chemotherapy or radiotherapy coupling.** Adjuvant therapy would shift
  the tumour axis. Irradiated beds are represented only through low perfusion.
- The outcome thresholds (80% bone, detection density 0.05, structural 50% at
  12 weeks, H₂ at 5× tolerance) are conventions and can be changed in `model.py`.
