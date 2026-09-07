---
output:
  pdf_document: default
  html_document: default
---
# Measurement Invariance Power App

Shiny app and R functions for computing analytic power of chi-square
difference tests of **scale-level** (parameter-total) vs **item-level**
measurement invariance in multiple-group one-factor CFA, generalizing the
code behind the manuscript in `Old Jordan Files for Claude/`.

## Run the app

```r
shiny::runApp("path/to/ShinyApp")
```

Requires: `lavaan`, `shiny`, `bslib`, `ggplot2`, `DT`.

## Files

| File | Contents |
|---|---|
| `app.R` | Shiny UI/server only (no statistics live here) |
| `R/models.R` | Population specification (`default_params`, `apply_noninvariance`, `build_structured_params`, `structured_sweep`, `pop_moments`) and lavaan model generation for any number of items *p* and groups *G* (`make_item_model`, `make_scale_model`) |
| `R/power.R` | `mi_power()` / `mi_power_sweep()`: MacCallum-Browne-Cai (2006) analytic power for metric/scalar/strict tests, both approaches, returned as tidy data frames |
| `R/plots.R` | `plot_power()`: ggplot2 power curves. Solid = scale-level, dashed = item-level, both in one panel so they stay comparable. Overlaid conditions are shades within each approach's hue (blue = scale-level, orange = item-level, dark = lowest value) and each curve is labelled with its value at the right-hand tip, so the plot stays readable with many conditions |
| `validation/validate_manuscript.R` | Reproduces every population RMSEA annotated in the *earlier* figure set (see the figure-map section) and checks power-formula equivalence with the original pipeline. Run: `Rscript validation/validate_manuscript.R` (all 44 checks pass) |
| `validation/test_app.R` | End-to-end checks of the app itself: UI structure, every sidebar control, the contents of each tab, the plot, the downloads, and the error paths. Run: `Rscript validation/test_app.R` (85 checks, about 30 seconds) |
| `marker_check.R` | Demonstrates that the scale-level test depends on which item is the marker while the item-level test does not |

## Scripting the functions directly

```r
invisible(lapply(list.files("R", full.names = TRUE), source))
library(lavaan)

# Manuscript (v10) Fig 4: k = 1..7 non-invariant loadings, +.1 each
conds <- structured_sweep(
  list(p = 8, n_groups = 2, type = "loadings", n_ni = 3, magnitude = .1),
  sweep_var = "n_ni", values = 1:7)
res <- mi_power_sweep(conds, tests = "metric")
plot_power(res$curves)
res$details   # population RMSEA, test df, N needed for 80% power
```

## Reproducing the manuscript figures in the app

Figure numbers refer to **Manuscript Version 10**
(`../Revision Spring 2025/Manuscript Version 10.docx`); the figures
themselves are archived in `Figures for Paper/`. Figure 1 is a schematic
of regression lines and is not produced by the app. Figures 2-6 are the
five power figures, and each is a single app plot: the colored curves are
the values of the input under "Draw one curve per value of", solid is
scale-level (the manuscript's *composite-level*) and dashed is item-level.

Common settings for all five (the app defaults unless noted):
2 groups, p = 8 items, base loading .7, base residual variance .51,
non-invariant parameters: **Loadings** ticked and the other two boxes
clear; direction: **All positive** (except Fig 2). The
significance level is fixed at alpha = .05. All five show only the metric
facet, so under "Tests to plot" tick Metric and untick Scalar and Strict.
N per group is 100-1000 for Figs 2-3 and 50-1000 for Figs 4-6.

The manuscript places non-invariance on different item positions than the
app (which always uses the last k items); because all other items are
exchangeable this yields identical fit statistics and power.

| Fig | Non-invariance settings | One curve per value of | What it shows |
|---|---|---|---|
| 2 | count: 6; magnitude: per-item; direction: **Totals cancel out**; N 100-1000 | Magnitude; values: `0.1, 0.2, 0.3` | Scale-level flat at .05 (all three solid curves coincide); item-level power climbs with magnitude |
| 3 | count: 6; magnitude: per-item; direction: All positive; N 100-1000 | Magnitude; values: `0.1, 0.2, 0.3` | Both approaches gain power; item-level somewhat higher |
| 4 | magnitude: per-item `.1`; direction: All positive | Number of non-invariant items; range 1-7 | Scale-level power increases monotonically in k; item-level peaks at k = 4 |
| 5 | count: 6; magnitude: **Total across non-invariant items** = `.6`; direction: All positive | Concentration of the total on one item; values: `0, 0.25, 0.5, 0.75, 1` | Scale-level unchanged (total is fixed); item-level rises as the total concentrates |
| 6 | count: 6; magnitude: per-item `.07`; direction: All positive | Number of items (p); values: `8, 12, 16, 20` | Scale-level power falls with scale length; item-level rises |

Endpoint powers at N = 1000 per group, read off the archived figures and
reproduced by the settings above (scale-level / item-level):

| Fig | Curve | Power |
|---|---|---|
| 2 | per-item .1, .2, .3 | .050/1.000 for all three (at N = 100: .050/.310, .050/.918, .050/1.000) |
| 3 | per-item .1, .2, .3 | .598/.682, .987/.999, 1.000/1.000 |
| 4 | k = 1 ... 7 | .067/.449, .118/.714, .205/.817, .322/.839, .459/.803, .598/.682, .724/.411 |
| 5 | largest item 17% ... 100% | .598/.682, .598/.976, then .599/1.000 |
| 6 | p = 8, 12, 16, 20 | .348/.362, .186/.640, .127/.718, .099/.743 |

Notes:
- **Terminology.** The manuscript says *composite-level*; the app, this
  README and the code say *scale-level*. They are the same thing: equality
  of the parameter totals across groups.
- **Figure 6's caption is off.** Version 10's caption and body text say the
  per-item difference is `.1` and the total is `.6`, but the plotted figure
  is `.07` per item (total `.42`) on 6 items; `.1` per item would put the
  p = 8 scale-level curve at .598 rather than the plotted .348. The body
  text also says "the total number of **invariant** items constant at 6",
  where the caption and the Discussion both mean *non-invariant*.
- The manuscript curves show power of the test matching the manipulated
  parameter type only. Left ticked, the scalar and strict facets are flat
  at alpha, because only loadings are manipulated in Figs 2-6.
- `validation/validate_manuscript.R` predates this figure set: its 44
  checks assert population RMSEA values for the **earlier** figures
  (`Old Jordan Files for Claude/`, "figures with captions.pdf" from
  `figure markdown4.Rmd`), whose numbering is unrelated to Version 10's.

## Method notes / differences from the original scripts

- Noncentrality is computed from the population minimum discrepancy:
  `nonc = (N_total - G) * (F0_constrained - F0_baseline)`, which for G = 2
  is numerically identical to the manuscript's `n * df * rmsea^2` but is
  correct for any G and independent of lavaan's multi-group RMSEA
  convention. Baselines are one invariance level lower (configural ->
  metric -> scalar -> strict), fit with the same approach.
- Test df are taken from the fitted models. This silently fixes a bug in
  the original `power.plot.doer.weak`/`.strict`, where the critical value
  used `nrow(pop.cov1-1)` = *p* df instead of *p - 1* for the metric test
  (a small effect, e.g. power .797 vs .819 at n = 300 in the earlier set's
  Fig 1a).
- Scale-level models for G groups impose G - 1 total-equality constraints
  per parameter class (plus the first-intercept anchor for scalar).
- Non-invariance can affect the last m groups (1 <= m <= G - 1), each
  receiving the same differences; the remaining groups equal group 1.
- "Concentration" redistributes a fixed total difference over the
  non-invariant items without changing that total: the largest item's share
  is 1/k + concentration * (1 - 1/k), the rest is split evenly over the
  other k - 1 items. 0 is the even split and 1 puts everything on one item,
  so at concentration = 1 the k-item condition reproduces the single-item
  condition exactly. This is the gradient in Version 10's Figure 5
  (item-level power rises as the total is concentrated, scale-level power
  stays put) while holding the number of non-invariant items fixed.
- Population latent variance is 1 and latent mean 0 in all groups;
  ad-hoc effect sizes (beta_spurious, d_spurious) are intentionally omitted.
