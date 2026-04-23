# TAAR1 Agonists — Triangulate Package Demo

**Author:** Alexandra Bannach-Brown  
**Project:** GALENOS living systematic review — animal-human triangulation workflow

---

## What this project does

This repository tests the [`triangulate`](https://mrcieu.r-universe.dev/triangulate) R package using a real GALENOS use-case: triangulating animal and human evidence on **trace amine-associated receptor 1 (TAAR1) agonists** (ulotaront, ralmitaront) for psychosis/schizophrenia.

The workflow takes study-level effect estimates from two evidence streams, attaches structured risk-of-bias and indirectness judgements, applies empirically anchored bias priors, and produces a combined forest plot showing original and bias-adjusted SMDs side-by-side with a risk-of-bias traffic-light grid.

---

## Clinical and animal context

- **Human RCTs:** Ulotaront and ralmitaront showed little difference versus placebo for overall symptoms in adults with acute schizophrenia. Four studies, 1291 participants; pooled SMD ~0.15 (95% CI −0.05 to 0.34).
- **Animal studies:** 188 effect sizes across 17 drug-level estimates from a multivariate meta-analysis (SYRCLE risk-of-bias framework). Pooled animal SMD ~1.11 (95% CI 0.84 to 1.38), reflecting reduced pro-psychotic drug-induced locomotor activity.

---

## Key source documents

| Resource | Link |
|---|---|
| `triangulate` package (Shapland et al 2024) | https://mrcieu.r-universe.dev/triangulate |
| `triangulate` preprint | https://www.medrxiv.org/content/10.1101/2024.09.20.24314046v1.full.pdf |
| TAAR1 living systematic review (Siafis et al 2024) | https://pmc.ncbi.nlm.nih.gov/articles/PMC11258611/ |
| GALENOS TAAR1 triangulation paper (Smith et al 2025) | https://www.cambridge.org/core/journals/the-british-journal-of-psychiatry/article/triangulating-evidence-from-the-galenos-living-systematic-review-on-trace-amineassociated-receptor-1-taar1-agonists-in-psychosis/B5202B48429F1E700F5274724FB176C0 |
| GALENOS GATE approach paper (Smith et al 2025) | https://www.cambridge.org/core/journals/the-british-journal-of-psychiatry/article/galenos-approach-to-triangulating-evidence-gate-transforming-the-landscape-of-psychiatric-research/94A657ABCCF4755661502785529339E0 |

---

## Repository structure

```
TAAR1agonist_triangulatePackage_demo/
├── TAAR1_triangulate_workflow.rmd   # Main analysis document (start here)
├── TAAR1_triangulate_workflow.html  # Rendered output
├── taar1-agonists.R                 # Full analysis script (sourced by the .rmd)
├── animal_df_full.csv               # 188-effect animal dataset (from GALENOS LSR3)
├── human_taar1.csv                  # Trial-level human RCT SMDs
├── taar1_drug_merged.csv            # 17 drug-level animal rows for triangulation plot
├── TAAR1_triangulate_plot.pdf       # Exported triangulation forest plot
└── Functions/
    ├── BiasCorrect_robviz.R         # Custom amended plotting function (main local adaptation)
    ├── helpers.R                    # ROB judgement/bias-type/direction cleaning
    ├── helpers-metafor.R            # p-value formatting and polygon annotation helpers
    ├── animal-helpers.R             # Multivariate animal model and forest-plot helpers
    └── rob_direction_modi.R         # Supporting plotting utilities
```

### How the .rmd and the .R script relate

`TAAR1_triangulate_workflow.rmd` does **not** duplicate the analysis code. It reads `taar1-agonists.R` with `readLines()` and injects code chunks directly from named markers in that file. This keeps a single source of truth: editing `taar1-agonists.R` automatically updates the rendered document on the next knit.

---

## Analysis workflow

### 1. Animal data — multivariate model (188 effects)
The visible triangulation plot shows 17 drug-level rows, but the animal summary polygon uses the original 188-effect multivariate model (`metafor::rma.mv()` with nested random effects). This model is recreated from `animal_df_full.csv` using `Functions/animal-helpers.R`.

### 2. Prepare data for `triangulate`
Both animal and human data are formatted so each row has:
- `yi` / `vi`: effect estimate and variance (SMD scale)
- `d1j … d10j`: risk-of-bias judgement per domain
- `d1t … d10t`: bias type (proportional / additive / none)
- `d1d … d10d`: bias direction (Away from null / Towards null / Unpredictable / None)

Animal data use **10 SYRCLE domains (D1–D10)**; human RCT data use **ROB2 domains (D1–D5)** with empty D6–D10 columns so both can be plotted on a shared grid.

### 3. `triangulate` package functions used
| Function | Purpose |
|---|---|
| `tri_to_long()` | Reshape wide domain columns to one row per result-domain pair |
| `tri_dat_check()` | Validate required columns |
| `tri_absolute_direction()` | Convert relative directions to absolute |
| `tri_append_bias()` | Join custom bias prior parameters |
| `tri_append_indirect()` | Join indirectness parameters (explicit "None" rows used here) |
| `tri_prep_data()` | Calculate `yi_adj` and `vi_adj` |

### 4. Custom bias priors
Priors are empirically anchored to each evidence stream's own meta-analytic uncertainty rather than package defaults.

**Animal priors** — anchored to the animal pooled 95% CI half-width (0.27 SMD):

| Judgement | Additive mean | Proportional mean |
|---|---:|---:|
| Low | 0.068 | 0.06 |
| Moderate | 0.135 | 0.12 |
| High | 0.27 | 0.24 |
| Critical | 0.54 | 0.48 |

**Human RCT priors** — anchored to the human pooled 95% CI half-width (0.195 SMD); proportional prior capped at 0.50 because the pooled effect is close to null:

| Judgement | Additive mean | Proportional mean |
|---|---:|---:|
| Low | 0.049 | 0.125 |
| Moderate | 0.098 | 0.25 |
| High | 0.195 | 0.50 |
| Critical | 0.39 | 1.00 |

### 5. Subgroup summaries
- **Animal:** `metafor::rma.mv()` fitted to the full 188-effect dataset, for both original and bias-adjusted panels.
- **Human:** `meta::metagen()` fitted to the trial-level rows, for both panels.
- Both are passed as **external summaries** to the plotting function (not fitted internally by the package).

### 6. Final plot
`Functions/BiasCorrect_robviz.R` contains a locally amended version of the package's internal `rob_direction()` plotting function. Key differences from the package default:

- Displays **D1–D10** (needed for SYRCLE animal domains)
- Keeps estimates on the **SMD (identity) scale** — no `exp()` transformation
- Separates confidence interval columns from the ROB traffic-light grid
- Accepts external subgroup summaries for original and adjusted panels separately
- Adds direction icons inside ROB squares (`<`/`>` for proportional, `v`/`^` for additive, `?` for unpredictable)
- Uses ASCII-safe text for reliable PDF output

---

## How to run

1. Open `TAAR1agonist_triangulatePackage_demo.Rproj` in RStudio.
2. Install required packages if needed:
   ```r
   install.packages(c("triangulate", "dplyr", "readr", "metafor", "meta", "ggplot2", "stringr"),
                    repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
   ```
3. Knit `TAAR1_triangulate_workflow.rmd` to reproduce `TAAR1_triangulate_workflow.html` and the PDF plot.

To run the script directly without knitting:
```r
source("taar1-agonists.R")
```

---

## Outputs

| File | Description |
|---|---|
| `TAAR1_triangulate_workflow.html` | Full rendered report with all code, tables, and the triangulation plot |
| `TAAR1_triangulate_plot.pdf` | Standalone triangulation forest plot (18 × 9 in) |

---

## Data provenance

- `animal_df_full.csv` — from the GALENOS LSR3 animal analysis: https://github.com/galenos-project/LSR3_taar1_A
- `human_taar1.csv` — trial-level SMDs extracted from the GALENOS TAAR1 living systematic review (Siafis et al 2024)
- `taar1_drug_merged.csv` — drug-level aggregation of the animal data used as visible rows in the triangulation plot
