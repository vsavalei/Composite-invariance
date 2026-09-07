# ---------------------------------------------------------------------------
# test_app.R -- end-to-end checks of the Shiny app: the UI structure, the
# server logic behind every sidebar control, the contents of each tab, the
# plot, the downloads, and the error paths.
#
# This file checks that the app wires the functions in R/ up correctly and
# displays the right thing; it does not re-derive the statistics themselves.
#
# Run from the ShinyApp folder:  Rscript validation/test_app.R
# ---------------------------------------------------------------------------

suppressMessages({
  library(shiny); library(bslib); library(DT); library(ggplot2); library(lavaan)
})
invisible(lapply(list.files("R", full.names = TRUE), source))

## ---- tiny test harness -----------------------------------------------------
n_pass <- 0L; n_fail <- 0L; failures <- character(0)

ok <- function(label, expr) {
  res <- tryCatch(isTRUE(expr), error = function(e) structure(FALSE, msg = conditionMessage(e)))
  if (isTRUE(res)) {
    n_pass <<- n_pass + 1L
    cat(sprintf("  ok    %s\n", label))
  } else {
    n_fail <<- n_fail + 1L
    msg <- attr(res, "msg")
    failures <<- c(failures, label)
    cat(sprintf("  FAIL  %s%s\n", label,
                if (is.null(msg)) "" else paste0("  [", msg, "]")))
  }
  invisible(NULL)
}
section <- function(x) cat(sprintf("\n== %s ==\n", x))

# rendered UI comes back in several pieces (tag html + dependencies)
ui_text <- function(x) paste(as.character(x), collapse = "")
# a population fitted at n = 1e7 is zero only to within numerical precision
ZERO <- 1e-6

## ---------------------------------------------------------------------------
## 1. UI structure
## ---------------------------------------------------------------------------
section("UI structure")

app_env <- new.env()
suppressMessages(sys.source("app.R", app_env))
tg <- htmltools::renderTags(app_env$ui)
html <- as.character(tg$html)
head_html <- paste(as.character(tg$head), collapse = "")
has_id <- function(id) grepl(paste0('id="', id, '"'), html, fixed = TRUE)

sidebar_inputs <- c("n_groups", "p", "base_loading", "base_residual",
                    "ni_loadings", "ni_intercepts", "ni_residuals",
                    "n_ni", "n_ni_groups", "mag_mode", "magnitude",
                    "magnitude_int", "magnitude_res", "pattern", "conc",
                    "compare_var",
                    "compare_k", "compare_g", "compare_p", "compare_mag",
                    "compare_conc", "include_metric", "include_scalar",
                    "include_strict", "n_range", "go")
ok("every sidebar input is rendered", all(vapply(sidebar_inputs, has_id, logical(1))))

outputs <- c("power_plot_ui", "details_scale", "details_item", "tbl_conditions",
             "tbl_loadings", "tbl_intercepts", "tbl_residuals",
             "dl_csv", "dl_png")
ok("every output/button is rendered", all(vapply(outputs, has_id, logical(1))))

ok("four tabs present",
   all(vapply(c("Power curves", "Details", "Population", "About"),
              function(x) grepl(paste0('data-value="', x, '"'), html, fixed = TRUE),
              logical(1))))

# a fillable card body turns the tables into flex items that overlap headings
bodies <- regmatches(html, gregexpr('<div class="card-body[^"]*"', html))[[1]]
ok("tab bodies are not fillable (no overlapping tables)",
   length(bodies) == 4 && !any(grepl("html-fill-container", bodies)))

narrow <- regmatches(html, gregexpr(
  '<div class="narrow-table">\\s*<div[^>]*id="(tbl_[a-z]+)"', html, perl = TRUE))[[1]]
ok("parameter tables are width-limited",
   setequal(sub('.*id="(tbl_[a-z]+)".*', "\\1", narrow),
            c("tbl_loadings", "tbl_intercepts", "tbl_residuals")) &&
   grepl(".narrow-table { max-width", head_html, fixed = TRUE))

# conditionalPanel expressions must reference inputs that exist
panel_conds <- unique(unlist(regmatches(
  html, gregexpr('data-display-if="[^"]*"', html))))
cond_inputs <- unique(unlist(regmatches(panel_conds,
  gregexpr("input\\.[A-Za-z_]+", panel_conds))))
cond_inputs <- sub("^input\\.", "", cond_inputs)
ok("conditional panels reference existing inputs",
   all(cond_inputs %in% sidebar_inputs))

ok("alpha is stated in the About tab, not an input",
   !has_id("alpha") && grepl("alpha = .05", html, fixed = TRUE))

## ---------------------------------------------------------------------------
## 2. Server: one testServer session driving every control
## ---------------------------------------------------------------------------
section("Server logic and tab contents")

defaults <- list(
  p = 8, n_groups = 2, base_loading = .7, base_residual = .51,
  ni_loadings = TRUE, ni_intercepts = FALSE, ni_residuals = FALSE,
  n_ni = 3, n_ni_groups = 1, mag_mode = "per_item",
  magnitude = .1, magnitude_int = .1, magnitude_res = .1,
  pattern = "uniform", conc = 0,
  compare_var = "none", compare_k = c(2, 5), compare_g = c(1, 1),
  compare_p = "4, 6, 8", compare_mag = "0.05, 0.1",
  compare_conc = "0, 0.5, 1", include_metric = TRUE, include_scalar = TRUE,
  include_strict = FALSE, n_range = c(50, 300))

testServer(shinyAppFile("app.R"), {
  go_n <- 0L
  run <- function(...) {                       # set inputs, then compute
    go_n <<- go_n + 1L
    do.call(session$setInputs, c(modifyList(defaults, list(...)), list(go = go_n)))
    results()
  }

  ## -- defaults --------------------------------------------------------------
  r <- run()
  ok("default run returns curves and details",
     !is.null(r) && nrow(r$curves) > 0 && nrow(r$details) == 4)
  ok("default run covers both approaches and both tests",
     setequal(r$curves$approach, c("scale", "item")) &&
     setequal(r$curves$test, c("metric", "scalar")))
  ok("N grid honours the slider and steps by 10",
     min(r$curves$n) == 50 && max(r$curves$n) == 300 &&
     all(diff(sort(unique(r$curves$n))) == 10))
  ok("power is a probability",
     all(r$curves$power >= 0 & r$curves$power <= 1))
  ok("details carry the four RMSEA columns",
     all(c("rmsea_configural", "rmsea_baseline", "rmsea", "rmsea_d") %in%
         names(r$details)))
  ok("configural model fits the population exactly",
     all(abs(r$details$rmsea_configural) < 1e-8))

  ## -- alpha is fixed at .05 -------------------------------------------------
  r <- run(n_ni = 3, magnitude = 0)
  ok("invariant population gives power = alpha = .05 everywhere",
     max(abs(r$curves$power - .05)) < 1e-8)

  ## -- non-invariant parameter classes: every combination --------------------
  # The three checkboxes are independent, so walk all eight subsets: each must
  # move exactly the ticked classes, and the population misfit of each test
  # must be non-zero for exactly the expected tests.
  combos <- expand.grid(loadings = c(FALSE, TRUE), intercepts = c(FALSE, TRUE),
                        residuals = c(FALSE, TRUE))
  mags <- c(loadings = .10, intercepts = .15, residuals = .05)
  cls_row <- c(loadings = "Loading", intercepts = "Intercept",
               residuals = "Residual variance")
  moved_ok <- detect_ok <- table_ok <- logical(nrow(combos))
  for (i in seq_len(nrow(combos))) {
    on <- names(combos)[unlist(combos[i, ])]
    r <- run(ni_loadings = combos$loadings[i], ni_intercepts = combos$intercepts[i],
             ni_residuals = combos$residuals[i], n_ni = 3,
             magnitude = mags[["loadings"]], magnitude_int = mags[["intercepts"]],
             magnitude_res = mags[["residuals"]],
             include_metric = TRUE, include_scalar = TRUE, include_strict = TRUE)
    prm <- r$params[[1]]
    # every class moves by its own magnitude on the last 3 items, or not at all
    moved_ok[i] <- all(vapply(names(mags), function(cl) {
      d <- prm[[cl]][, 2] - prm[[cl]][, 1]
      want <- if (cl %in% on) mags[[cl]] else 0
      all(abs(d[6:8] - want) < 1e-12) && all(d[1:5] == 0)
    }, logical(1)))
    # The metric and scalar tests have non-zero population misfit only for
    # their own parameter class. The strict test also has non-zero misfit
    # whenever the loadings are non-invariant, even when the population
    # residual variances are equal across groups: the metric and scalar
    # models constrain the loadings to be equal, and because the population
    # loadings are not equal, the estimated residual variances differ across
    # groups in those models. Constraining them equal at the strict level
    # therefore adds misfit. Non-invariant intercepts do not have this
    # effect, because intercepts enter only the mean structure while the
    # residual variances are estimated from the covariance structure.
    misfit <- c(metric = "loadings" %in% on,
                scalar = "intercepts" %in% on,
                strict = any(c("loadings", "residuals") %in% on))
    detect_ok[i] <- all(vapply(names(misfit), function(tst) {
      rd <- r$details$rmsea_d[r$details$test == tst]
      length(rd) == 2 && if (misfit[[tst]]) all(rd > ZERO) else all(rd < ZERO)
    }, logical(1)))
    cond <- conditions_data()
    table_ok[i] <- if (!length(on)) is.null(cond) else
      setequal(cond$Parameter, unname(cls_row[on]))
  }
  ok("every combination moves exactly the ticked parameter classes",
     all(moved_ok))
  ok("every combination is detected by exactly the matching tests",
     all(detect_ok))
  ok("the conditions table lists exactly the ticked parameter classes",
     all(table_ok))

  r <- run(ni_loadings = FALSE, ni_intercepts = FALSE, ni_residuals = FALSE,
           include_strict = TRUE)
  ok("no class ticked leaves every test at alpha = .05",
     !is.null(r) && max(abs(r$curves$power - .05)) < 1e-8)

  # the one effect that crosses parameter classes, asserted separately so
  # that a change to it fails a named check
  r <- run(ni_loadings = TRUE, ni_intercepts = FALSE, ni_residuals = FALSE,
           include_strict = TRUE)
  ok("non-invariant loadings give the strict test non-zero misfit",
     all(r$details$rmsea_d[r$details$test == "strict"] > ZERO))
  r <- run(ni_loadings = FALSE, ni_intercepts = TRUE, ni_residuals = FALSE,
           include_strict = TRUE)
  ok("non-invariant intercepts give the strict test zero misfit",
     all(r$details$rmsea_d[r$details$test == "strict"] < ZERO))

  # combinations that the old single-choice control could not express
  r <- run(ni_loadings = TRUE, ni_intercepts = FALSE, ni_residuals = TRUE,
           magnitude = .1, magnitude_res = .05, include_strict = TRUE)
  ok("loadings + residual variances: metric and strict detect, scalar does not",
     all(r$details$rmsea_d[r$details$test %in% c("metric", "strict")] > ZERO) &&
     all(r$details$rmsea_d[r$details$test == "scalar"] < ZERO))
  r <- run(ni_loadings = FALSE, ni_intercepts = TRUE, ni_residuals = TRUE,
           magnitude_int = .15, magnitude_res = .05, include_strict = TRUE)
  ok("intercepts + residual variances: scalar and strict detect, metric does not",
     all(r$details$rmsea_d[r$details$test %in% c("scalar", "strict")] > ZERO) &&
     all(r$details$rmsea_d[r$details$test == "metric"] < ZERO))

  ok("unticking a class stops it being applied",
     {
       a <- run(ni_loadings = TRUE, ni_residuals = TRUE, magnitude_res = .05)
       b <- run(ni_loadings = TRUE, ni_residuals = FALSE, magnitude_res = .05)
       pa <- a$params[[1]]; pb <- b$params[[1]]
       identical(pa$loadings, pb$loadings) &&
         any(pa$residuals[, 2] != pa$residuals[, 1]) &&
         all(pb$residuals[, 2] == pb$residuals[, 1])
     })

  # a magnitude only matters for a class that is ticked
  ok("a magnitude for an unticked class is ignored",
     identical(run(ni_loadings = TRUE, ni_intercepts = FALSE,
                   magnitude_int = .9)$params[[1]],
               run(ni_loadings = TRUE, ni_intercepts = FALSE,
                   magnitude_int = .15)$params[[1]]))

  r <- run(ni_loadings = TRUE, ni_intercepts = TRUE, ni_residuals = TRUE,
           magnitude_res = -.9, include_strict = TRUE)
  ok("a residual magnitude that drives variances negative is rejected",
     is.null(r))

  # a magnitude sweep moves every ticked class together
  r <- run(ni_loadings = TRUE, ni_intercepts = TRUE, ni_residuals = TRUE,
           compare_var = "magnitude", compare_mag = "0.05, 0.2",
           include_strict = TRUE)
  ok("a magnitude sweep moves every ticked class to the swept value",
     length(r$params) == 2 &&
     all(vapply(names(r$params), function(nm) {
       v <- as.numeric(sub(".*= ", "", nm))
       prm <- r$params[[nm]]
       all(vapply(c("loadings", "intercepts", "residuals"), function(cl)
         abs(prm[[cl]][8, 2] - prm[[cl]][8, 1] - v) < 1e-12, logical(1)))
     }, logical(1))))

  ## -- direction patterns ----------------------------------------------------
  r <- run(pattern = "uniform", n_ni = 4)
  ok("all positive: scale-level test has a violation",
     all(r$details$rmsea_d[r$details$test == "metric" &
                           r$details$approach == "scale"] > ZERO))
  for (pat in c("cancel", "cancel_one")) {
    r <- run(pattern = pat, n_ni = 5)
    d <- r$details[r$details$test == "metric", ]
    ok(sprintf("%s: totals cancel (scale-level RMSEA_d is zero)", pat),
       d$rmsea_d[d$approach == "scale"] < ZERO)
    ok(sprintf("%s: item-level still detects the differences", pat),
       d$rmsea_d[d$approach == "item"] > ZERO)
    ok(sprintf("%s: power stays at alpha for the scale-level test", pat),
       max(abs(r$curves$power[r$curves$test == "metric" &
                              r$curves$approach == "scale"] - .05)) < 1e-4)
    ok(sprintf("%s: the population totals are exactly equal", pat),
       abs(diff(colSums(r$params[[1]]$loadings))) < 1e-12)
  }

  ## -- concentration ---------------------------------------------------------
  r <- run(compare_var = "concentration", n_ni = 6, compare_conc = "0, 0.5, 1")
  d <- r$details[r$details$test == "metric", ]
  sc <- d[d$approach == "scale", ]; it <- d[d$approach == "item", ]
  ok("concentration: three curves", length(unique(r$curves$condition)) == 3)
  ok("concentration: total is held fixed, so scale-level is unchanged",
     diff(range(sc$rmsea_d)) < 1e-3)
  ok("concentration: item-level power rises with concentration",
     all(diff(it$rmsea_d) > 0))
  ok("concentration: curves labelled by the largest item's share",
     all(grepl("largest item = .*% of total", unique(r$curves$condition))))

  ## -- comparison variables --------------------------------------------------
  r <- run(compare_var = "n_ni", compare_k = c(2, 5))
  ok("compare items: one curve per value in the range",
     setequal(unique(r$curves$condition),
              paste(2:5, "non-invariant")))
  r <- run(compare_var = "p", compare_p = "6, 8, 10")
  ok("compare p: one curve per listed value",
     setequal(unique(r$curves$condition), paste("p =", c(6, 8, 10))))
  r <- run(compare_var = "magnitude", compare_mag = "0.05, 0.2")
  ok("compare magnitude: labelled per item",
     setequal(unique(r$curves$condition),
              c("per-item diff = 0.05", "per-item diff = 0.2")))
  r <- run(compare_var = "magnitude", mag_mode = "total", compare_mag = "0.3, 0.6")
  ok("compare magnitude: labelled as a total when the mode says so",
     all(grepl("^total diff", unique(r$curves$condition))))
  r <- run(n_groups = 4, n_ni_groups = 1, compare_var = "n_ni_groups",
           compare_g = c(1, 3))
  ok("compare non-invariant groups: one curve per value",
     setequal(unique(r$curves$condition),
              paste("non-invariant groups =", 1:3)))

  ## -- groups ----------------------------------------------------------------
  for (G in 2:4) {
    r <- run(n_groups = G, n_ni_groups = 1, n_ni = 3)
    d <- r$details[r$details$test == "metric", ]
    ok(sprintf("G = %d: scale-level test has G - 1 df", G),
       d$df_test[d$approach == "scale"] == G - 1)
    ok(sprintf("G = %d: item-level test has (p - 1)(G - 1) df", G),
       d$df_test[d$approach == "item"] == 7 * (G - 1))
  }

  ## -- test selection --------------------------------------------------------
  r <- run(include_metric = TRUE, include_scalar = FALSE, include_strict = FALSE)
  ok("metric only", setequal(r$curves$test, "metric"))
  r <- run(include_metric = FALSE, include_scalar = TRUE, include_strict = TRUE)
  ok("scalar and strict only", setequal(r$curves$test, c("scalar", "strict")))
  r <- run(include_metric = TRUE, include_scalar = TRUE, include_strict = TRUE)
  ok("all three tests", setequal(r$curves$test, c("metric", "scalar", "strict")))
  ok("plot area grows for three tests",
     grepl("840px", ui_text(output$power_plot_ui), fixed = TRUE))
  r <- run(include_metric = TRUE, include_scalar = TRUE)
  ok("plot area is shorter for two tests",
     grepl("470px", ui_text(output$power_plot_ui), fixed = TRUE))

  ## -- Details tab -----------------------------------------------------------
  # DT ships its rows over Ajax, so the rendered widget only carries the
  # headers: check the displayed values through the functions that build them
  r <- run(compare_var = "none", n_ni = 3, magnitude = .1, n_range = c(50, 1000))
  ds <- details_data("scale"); di <- details_data("item")
  ok("Details: both tables render",
     nchar(ui_text(output$details_scale)) > 0 &&
     nchar(ui_text(output$details_item)) > 0)
  ok("Details: expected columns",
     identical(names(di), c("Test", "RMSEA (configural)", "RMSEA (model)",
                            "RMSEA_d", "df (test)", "N for 60% power",
                            "N for 80% power")))
  ok("Details: one row per test", nrow(ds) == 2 && nrow(di) == 2)
  ok("Details: unreached 80% power shown as '> max', not blank",
     all(ds$`N for 80% power` == "> 1000"))
  ok("Details: a reached 80% power shows the sample size",
     di$`N for 80% power`[di$Test == "Metric"] == "970")
  ok("Details: unreached 60% power shown as '> max', not blank",
     all(ds$`N for 60% power` == "> 1000"))
  ok("Details: 60% power is reached at a smaller N than 80% power",
     as.numeric(di$`N for 60% power`[di$Test == "Metric"]) <
       as.numeric(di$`N for 80% power`[di$Test == "Metric"]))
  ok("Details: configural RMSEA column is zero for this population",
     all(ds$`RMSEA (configural)` == 0))
  ok("Details: no Curve column for a single setting",
     !("Curve" %in% names(ds)))
  r <- run(compare_var = "n_ni", compare_k = c(2, 4))
  ok("Details: Curve column appears when curves are compared",
     "Curve" %in% names(details_data("scale")) &&
     setequal(details_data("scale")$Curve, paste(2:4, "non-invariant")))

  ## -- Population tab --------------------------------------------------------
  r <- run(compare_var = "concentration", n_ni = 6, magnitude = .1,
           compare_conc = "0, 1")
  pl <- r$params
  ok("Population: the run stores one population per curve",
     !is.null(pl) && length(pl) == 2)
  tot <- vapply(pl, function(prm) sum(prm$loadings[, 2] - prm$loadings[, 1]), 0)
  ok("Population: concentration keeps the total difference at .6",
     all(abs(tot - .6) < 1e-12))
  cond <- conditions_data()
  ok("Population: conditions table renders",
     nchar(ui_text(output$tbl_conditions)) > 0)
  ok("Population: one row per curve, with items and a Total column",
     nrow(cond) == 2 && all(c("Parameter", "Curve", "y1", "y8", "Total") %in%
                            names(cond)))
  ok("Population: the listed differences add up to the stated total",
     all(abs(cond$Total - .6) < 1e-8))
  ok("Population: differences are the ones the settings imply",
     all(cond[1, paste0("y", 3:8)] == .1) &&
     cond[2, "y8"] == .6 && all(cond[2, paste0("y", 3:7)] == 0))
  ok("Population: only parameters that differ are listed",
     all(cond$Parameter == "Loading"))
  ok("Population: Curve column names the compared curves",
     setequal(cond$Curve, c("largest item = 17% of total",
                            "largest item = 100% of total")))
  ok("Population: parameter tables follow the sidebar",
     sum(population()$loadings[, 2] - population()$loadings[, 1] > 0) == 6)
  ok("Population: parameter tables render",
     nchar(ui_text(output$tbl_loadings)) > 0 &&
     nchar(ui_text(output$tbl_intercepts)) > 0 &&
     nchar(ui_text(output$tbl_residuals)) > 0)
  run(compare_var = "none", ni_loadings = TRUE, n_ni = 2, magnitude = .3)
  ok("Population: the tables track a later sidebar change",
     sum(population()$loadings[, 2] - population()$loadings[, 1] > 0) == 2 &&
     abs(population()$loadings[8, 2] - population()$loadings[8, 1] - .3) < 1e-12)
  ok("Population: the tables carry the item and group names",
     identical(rownames(population()$loadings), paste0("y", 1:8)) &&
     identical(colnames(population()$loadings), paste0("group", 1:2)))

  ## -- error paths -----------------------------------------------------------
  r <- run(include_metric = FALSE, include_scalar = FALSE,
           include_strict = FALSE)
  ok("no test selected clears the result", is.null(r))
  r <- run(compare_var = "p", compare_p = "not a number")
  ok("unparseable comparison list clears the result", is.null(r))
  r <- run(compare_var = "n_ni", compare_k = c(1, 5), pattern = "cancel")
  ok("impossible values are skipped, the rest still run",
     length(unique(r$curves$condition)) == 4 &&
     !("1 non-invariant" %in% r$curves$condition))
  r <- run(compare_var = "n_ni_groups", n_groups = 2, compare_g = c(1, 1))
  ok("comparing groups with G = 2 still draws the single possible curve",
     !is.null(r) && length(unique(r$curves$condition)) == 1)

  ## -- downloads -------------------------------------------------------------
  r <- run(compare_var = "none", n_ni = 3, n_range = c(50, 200))
  ok("Population: no Curve column for a single setting",
     !("Curve" %in% names(conditions_data())))
  csv <- tempfile(fileext = ".csv")
  write.csv(r$curves, csv, row.names = FALSE)
  back <- read.csv(csv)
  ok("CSV download content round-trips",
     nrow(back) == nrow(r$curves) &&
     all(c("n", "power", "test", "approach", "condition") %in% names(back)))
  unlink(csv)
})

## ---------------------------------------------------------------------------
## 2b. Non-invariant parameter classes, at the function level
## ---------------------------------------------------------------------------
section("Non-invariant parameter classes")

ok("the three classes are named in a fixed order",
   identical(MI_PARAM_CLASSES, c("loadings", "intercepts", "residuals")))
ok("a selection is put in the canonical order",
   identical(mi_param_classes(c("residuals", "loadings")),
             c("loadings", "residuals")))
ok("repeats collapse", identical(mi_param_classes(c("loadings", "loadings")),
                                 "loadings"))
ok("an empty selection is allowed",
   identical(mi_param_classes(character(0)), character(0)) &&
   identical(mi_param_classes(NULL), character(0)))
ok("the old 'both' shorthand still means loadings + intercepts",
   identical(mi_param_classes("both"), c("loadings", "intercepts")))
ok("the old 'all' shorthand still means every class",
   identical(mi_param_classes("all"), MI_PARAM_CLASSES))
ok("an unknown class is rejected",
   inherits(tryCatch(mi_param_classes("slopes"), error = function(e) e), "error"))

bsp <- function(type, n_ni = 3, ...)
  build_structured_params(p = 8, n_groups = 2, base_loading = .7,
                          base_residual = .51, type = type, n_ni = n_ni,
                          magnitude = .1, magnitude_intercept = .15,
                          magnitude_residual = .05, ...)
moved <- function(prm) {
  m <- vapply(MI_PARAM_CLASSES,
              function(cl) max(abs(prm[[cl]][, 2] - prm[[cl]][, 1])), 0)
  names(m)[m > 0]
}
subsets <- list(character(0), "loadings", "intercepts", "residuals",
                c("loadings", "intercepts"), c("loadings", "residuals"),
                c("intercepts", "residuals"), MI_PARAM_CLASSES)
ok("every subset moves exactly the classes it names",
   all(vapply(subsets, function(s) setequal(moved(bsp(s)), s), logical(1))))
ok("each class takes its own magnitude",
   {
     prm <- bsp(MI_PARAM_CLASSES)
     abs(prm$loadings[8, 2]   - prm$loadings[8, 1]   - .10) < 1e-12 &&
     abs(prm$intercepts[8, 2] - prm$intercepts[8, 1] - .15) < 1e-12 &&
     abs(prm$residuals[8, 2]  - prm$residuals[8, 1]  - .05) < 1e-12
   })
ok("a single class still uses its own magnitude, not the loading one",
   abs(bsp("residuals")$residuals[8, 2] -
       bsp("residuals")$residuals[8, 1] - .05) < 1e-12)
ok("order of the requested classes does not change the population",
   identical(bsp(c("residuals", "intercepts", "loadings")),
             bsp(MI_PARAM_CLASSES)))
ok("the shorthands build the same populations as the explicit subsets",
   identical(bsp("both"), bsp(c("loadings", "intercepts"))) &&
   identical(bsp("all"), bsp(MI_PARAM_CLASSES)))
ok("an empty selection is a fully invariant population",
   identical(bsp(character(0)),
             default_params(8, 2, .7, base_residual = .51)))
ok("n_ni = 0 is invariant whatever is selected",
   identical(bsp(MI_PARAM_CLASSES, n_ni = 0),
             default_params(8, 2, .7, base_residual = .51)))
ok("a total magnitude is split within each class separately",
   {
     prm <- bsp(MI_PARAM_CLASSES, magnitude_is_total = TRUE)
     abs(sum(prm$loadings[, 2]   - prm$loadings[, 1])   - .10) < 1e-12 &&
     abs(sum(prm$intercepts[, 2] - prm$intercepts[, 1]) - .15) < 1e-12 &&
     abs(sum(prm$residuals[, 2]  - prm$residuals[, 1])  - .05) < 1e-12
   })
ok("concentration holds each class's own total",
   {
     tot <- function(cc) {
       prm <- bsp(MI_PARAM_CLASSES, magnitude_is_total = TRUE, concentration = cc)
       vapply(MI_PARAM_CLASSES,
              function(cl) sum(prm[[cl]][, 2] - prm[[cl]][, 1]), 0)
     }
     max(abs(tot(0) - tot(1))) < 1e-12
   })
ok("a cancelling direction zeroes the total of every selected class",
   {
     prm <- bsp(MI_PARAM_CLASSES, n_ni = 4, pattern = "cancel")
     all(vapply(MI_PARAM_CLASSES,
                function(cl) abs(sum(prm[[cl]][, 2] - prm[[cl]][, 1])) < 1e-12,
                logical(1)))
   })
ok("a selection that would drive a residual variance negative is rejected",
   inherits(tryCatch(bsp(MI_PARAM_CLASSES, magnitude_residual = -.9),
                     error = function(e) e), "error"))
ok("non-invariance reaches the last n_ni_groups of several groups",
   {
     prm <- build_structured_params(p = 8, n_groups = 4, base_residual = .51,
                                    type = c("loadings", "residuals"),
                                    n_ni = 3, magnitude = .1,
                                    magnitude_residual = .05, n_ni_groups = 2)
     all(prm$loadings[, 1] == prm$loadings[, 2]) &&
       all(prm$loadings[, 3] == prm$loadings[, 4]) &&
       abs(prm$loadings[8, 3] - prm$loadings[8, 1] - .1) < 1e-12 &&
       abs(prm$residuals[8, 3] - prm$residuals[8, 1] - .05) < 1e-12 &&
       all(prm$intercepts[, 3] == prm$intercepts[, 1])
   })
ok("a sweep over a combination keeps the combination in every condition",
   {
     s <- structured_sweep(list(p = 8, n_groups = 2, base_residual = .51,
                                type = c("loadings", "residuals"), n_ni = 3,
                                magnitude = .1, magnitude_residual = .05),
                           "n_ni", 2:4)
     length(s) == 3 &&
       all(vapply(s, function(prm)
         setequal(moved(prm), c("loadings", "residuals")), logical(1)))
   })
ok("a magnitude sweep sets every class to the swept value",
   {
     s <- structured_sweep(list(p = 8, n_groups = 2, base_residual = .51,
                                type = MI_PARAM_CLASSES, n_ni = 3,
                                magnitude = .1, magnitude_intercept = .15,
                                magnitude_residual = .05),
                           "magnitude", c(.02, .2))
     all(vapply(seq_along(s), function(i) {
       v <- c(.02, .2)[i]
       all(vapply(MI_PARAM_CLASSES, function(cl)
         abs(s[[i]][[cl]][8, 2] - s[[i]][[cl]][8, 1] - v) < 1e-12, logical(1)))
     }, logical(1)))
   })

## ---------------------------------------------------------------------------
## 3. Plot
## ---------------------------------------------------------------------------
section("Power curve plot")

prm <- build_structured_params(8, 2, type = "loadings", n_ni = 3, magnitude = .1)
res1 <- mi_power(prm, tests = c("metric", "scalar"), n_seq = seq(100, 500, 10))
g1 <- plot_power(res1$curves, ncol = 2)
b1 <- ggplot2::ggplot_build(g1)

ok("plot builds for a single setting", inherits(g1, "ggplot"))
ok("one facet per test", length(unique(b1$data[[2]]$PANEL)) == 2)
ok("solid = scale-level, dashed = item-level",
   setequal(b1$data[[2]]$linetype, c("solid", "42")))
ok("y axis is the full 0-1 probability range",
   all(abs(b1$layout$panel_params[[1]]$y.range - c(0, 1)) < .02))
ok("x axis starts exactly at the smallest N",
   abs(min(b1$layout$panel_params[[1]]$x.range) - 100) < 1e-8)
ok("the smallest N is labelled",
   100 %in% b1$layout$panel_params[[1]]$x$breaks)
scale_name <- function(g, aes) {
  s <- g$scales$get_scales(aes)
  if (is.null(s) || inherits(s$name, "waiver")) NULL else s$name
}
ok("no legend title when approach is the only distinction",
   is.null(scale_name(g1, "colour")))

cl <- structured_sweep(list(p = 8, n_groups = 2, type = "loadings", n_ni = 3,
                            magnitude = .1), "magnitude", c(.05, .1, .2))
res2 <- mi_power_sweep(cl, tests = c("metric", "scalar"),
                       n_seq = seq(50, 300, 10))
g2 <- plot_power(res2$curves, ncol = 2, color_title = "Magnitude")
b2 <- ggplot2::ggplot_build(g2)
ok("plot builds for compared curves", inherits(g2, "ggplot"))
ok("one line per curve per approach, drawn in every facet",
   length(unique(b2$data[[2]]$group)) == 3 * 2 &&
   nrow(unique(b2$data[[2]][, c("group", "PANEL")])) == 3 * 2 * 2)
ok("no text layer when there is nothing to compare", length(b1$data) == 2)

# both approaches share one panel, so each keeps its own hue and the compared
# values are shades within it
lines2 <- b2$data[[2]]
ok("scale-level curves are shades of the scale-level hue",
   setequal(unique(lines2$colour[lines2$linetype == "solid"]),
            mi_condition_colors("scale", 3)))
ok("item-level curves are shades of the item-level hue",
   setequal(unique(lines2$colour[lines2$linetype == "42"]),
            mi_condition_colors("item", 3)))
ok("every compared value gets its own shade",
   length(unique(lines2$colour)) == 3 * 2)
ok("shades run dark to light with the compared value",
   identical(mi_condition_colors("scale", 3)[1], "#052F5F"))

tips2 <- b2$data[[3]]
ok("one tip label per curve in every facet", nrow(tips2) == 3 * 2 * 2)
ok("tip labels are trimmed to the compared value",
   setequal(unique(tips2$label), c("0.05", "0.1", "0.2")))
ok("tip labels stay inside the panel", all(tips2$y >= 0 & tips2$y <= 1))
ok("tip labels take the colour of their curve",
   setequal(unique(tips2$colour), unique(lines2$colour)))
ok("the x axis leaves room for the tip labels",
   diff(b2$layout$panel_params[[1]]$x.range) > diff(range(res2$curves$n)) * 1.05)

ok("colour scale still records what varies",
   identical(scale_name(g2, "colour"), "Magnitude"))
ok("the colour key is dropped in favour of the tip labels",
   identical(g2$scales$get_scales("colour")$guide, "none"))
ok("the caption names what the tip labels are",
   grepl("magnitude", g2$labels$caption, fixed = TRUE))
ok("linetype legend is titled for compared curves",
   identical(scale_name(g2, "linetype"), "Approach"))

## tip-label helpers ---------------------------------------------------------
ok("tip labels trim every sweep label format",
   identical(mi_tip_labels(c("p = 12", "3 non-invariant",
                             "largest item = 17% of total",
                             "per-item diff = 0.05",
                             "non-invariant groups = 2")),
             c("12", "3", "17%", "0.05", "2")))
ok("tip labels keep the full name when trimming would be ambiguous",
   identical(mi_tip_labels(c("p = 8 of 2", "p = 8 of 3")),
             c("p = 8 of 2", "p = 8 of 3")))
ok("tip labels keep names that carry no number",
   identical(mi_tip_labels(c("custom", "other")), c("custom", "other")))

# a long sweep crowds the tips; they must stay ordered and on the panel
cl_big <- structured_sweep(list(p = 8, n_groups = 2, type = "loadings",
                                n_ni = 6, magnitude = .1),
                           "p", seq(8, 30, 2))
res3 <- mi_power_sweep(cl_big, tests = "metric", n_seq = seq(50, 1000, 50))
b3 <- ggplot2::ggplot_build(plot_power(res3$curves, color_title = "Number of items"))
tips3 <- b3$data[[3]]
ok("crowded tip labels all stay on the panel",
   nrow(tips3) == 12 * 2 && all(tips3$y >= 0 & tips3$y <= 1))
ok("crowded tip labels keep the order of the curves they label",
   !is.unsorted(rev(tips3$y[order(-tips3$y)])) &&
   !any(duplicated(round(tips3$y, 6))))

png_file <- tempfile(fileext = ".png")
suppressMessages(ggsave(png_file, g2, width = 11, height = 5, dpi = 100))
ok("plot saves as PNG", file.exists(png_file) && file.size(png_file) > 5000)
unlink(png_file)

## ---------------------------------------------------------------------------
section("Summary")
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) {
  cat("Failed checks:\n"); cat(paste0("  - ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("ALL APP TESTS PASSED\n")
