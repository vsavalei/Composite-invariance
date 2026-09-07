# ---------------------------------------------------------------------------
# Shiny app: power of scale-level vs item-level measurement invariance tests
#
# All statistics live in R/models.R and R/power.R (sourced below); this file
# only collects inputs, dispatches to those functions, and renders output.
# ---------------------------------------------------------------------------

library(shiny)
library(bslib)
library(ggplot2)
library(DT)
suppressMessages(library(lavaan))

invisible(lapply(list.files("R", full.names = TRUE), source))

num_list <- function(txt) {
  v <- suppressWarnings(as.numeric(strsplit(txt, "[,; ]+")[[1]]))
  v[!is.na(v)]
}

TEST_NAMES <- c(metric = "Metric", scalar = "Scalar", strict = "Strict")

# legend heading naming what the compared curves vary in
COMPARE_TITLES <- c(n_ni = "Non-invariant items",
                    n_ni_groups = "Non-invariant groups",
                    p = "Number of items",
                    magnitude = "Magnitude",
                    concentration = "Share on one item")

ui <- page_sidebar(
  title = "Measurement Invariance Power: Scale-Level vs Item-Level Tests",
  theme = bs_theme(bootswatch = "flatly", font_scale = 0.75),
  tags$head(tags$style(HTML(
    # parameter tables hold short numbers: keep them compact instead of
    # stretching every column across the panel
    ".narrow-table { max-width: 560px; }
     .narrow-table table.dataTable { width: auto !important; }
     .narrow-table table.dataTable th,
     .narrow-table table.dataTable td { padding: 3px 10px; text-align: right; }
     table.dataTable th, table.dataTable td { white-space: nowrap; }"))),
  sidebar = sidebar(
    width = 330,
    sliderInput("n_groups", "Number of groups", 2, 6, 2, step = 1, ticks = FALSE),
    conditionalPanel(
      "input.compare_var != 'p'",
      sliderInput("p", "Number of items (p)", 4, 20, 8, step = 1, ticks = FALSE)),
    conditionalPanel(
      "input.compare_var == 'p'",
      helpText("Number of items: one curve per value entered under ",
               "'Draw one curve per value of' below.")),
    numericInput("base_loading", "Base factor loading", .7, .1, 2, .05),
    numericInput("base_residual", "Base residual variance", .51, .01, 5, .01),
    helpText("Latent variance is 1 and latent mean 0 in every population ",
             "group. Sample sizes are equal across groups."),
    tags$label("Non-invariant parameters (added to the last group(s))",
               class = "control-label"),
    checkboxInput("ni_loadings", "Loadings", TRUE),
    checkboxInput("ni_intercepts", "Intercepts", FALSE),
    checkboxInput("ni_residuals", "Residual variances", FALSE),
    conditionalPanel(
      "!input.ni_loadings && !input.ni_intercepts && !input.ni_residuals",
      helpText("Nothing is non-invariant, so every test will sit at its ",
               "Type I error rate of .05.")),
    conditionalPanel(
      "input.compare_var != 'n_ni'",
      sliderInput("n_ni", "Number of non-invariant items", 1, 7, 3,
                  step = 1, ticks = FALSE)),
    conditionalPanel(
      "input.compare_var == 'n_ni'",
      helpText("Number of non-invariant items: one curve per value in ",
               "the range chosen under 'Draw one curve per value of' below.")),
    conditionalPanel(
      "input.n_groups > 2 && input.compare_var != 'n_ni_groups'",
      sliderInput("n_ni_groups", "Number of non-invariant groups", 1, 1, 1,
                  step = 1, ticks = FALSE)),
    conditionalPanel(
      "input.compare_var == 'n_ni_groups'",
      helpText("Number of non-invariant groups: one curve per value in ",
               "the range chosen under 'Draw one curve per value of' below.")),
    radioButtons("mag_mode", "Magnitude specified as",
                 c("Per-item difference" = "per_item",
                   "Total across non-invariant items" = "total")),
    conditionalPanel(
      "input.compare_var != 'magnitude'",
      # loadings, intercepts and residual variances are on different scales,
      # so each selected class gets its own difference
      conditionalPanel(
        "input.ni_loadings",
        numericInput("magnitude", "Magnitude for loading differences", .1,
                     step = .01)),
      conditionalPanel(
        "input.ni_intercepts",
        numericInput("magnitude_int", "Magnitude for intercept differences", .1,
                     step = .01)),
      conditionalPanel(
        "input.ni_residuals",
        numericInput("magnitude_res",
                     "Magnitude for residual variance differences", .1,
                     step = .01)),
      conditionalPanel(
        "input.ni_loadings || input.ni_intercepts || input.ni_residuals",
        helpText("Each magnitude is added in every non-invariant group."))),
    conditionalPanel(
      "input.compare_var == 'magnitude'",
      helpText("Magnitude: one curve per value entered under ",
               "'Draw one curve per value of' below.")),
    radioButtons("pattern", "Direction of differences",
                 c("All positive" = "uniform",
                   "Totals cancel out (differences spread across the items in opposite directions)" = "cancel",
                   "All positive, one negative (totals cancel out)" = "cancel_one")),
    conditionalPanel(
      "input.pattern == 'uniform' && input.compare_var != 'concentration'",
      sliderInput("conc", "Concentration of the total on one item", 0, 1, 0,
                  step = .05, ticks = FALSE)),
    conditionalPanel(
      "input.pattern == 'uniform' && input.compare_var == 'concentration'",
      helpText("Concentration: one curve per value entered under ",
               "'Draw one curve per value of' below.")),
    conditionalPanel(
      "input.pattern == 'uniform'",
      helpText("Concentration holds the total difference fixed and changes ",
               "how it is spread: 0 splits it evenly over the non-invariant ",
               "items, 1 puts all of it on a single item and leaves the ",
               "others invariant.")),
    conditionalPanel(
      "input.pattern != 'uniform'",
      helpText("Cancelling directions need at least 2 non-invariant items. ",
               "'Spread' uses alternating signs, recentered so the totals ",
               "are exactly equal (with an even count this is simply ",
               "+/- the magnitude); 'one negative' gives the last item ",
               "minus (count - 1) times the magnitude.")),
    helpText("Differences are added to the last k items of each ",
             "non-invariant group (always the last group(s), all with the ",
             "same differences); all other groups equal the first group. ",
             "Item 1 is always invariant: its loading is fixed to 1 in ",
             "every group in the analysis models, and in scalar and strict ",
             "models its intercept is also constrained equal across groups."),
    selectInput("compare_var", "Draw one curve per value of",
                c("None (single setting)" = "none",
                  "Number of non-invariant items" = "n_ni",
                  "Number of non-invariant groups" = "n_ni_groups",
                  "Number of items (p)" = "p",
                  "Magnitude" = "magnitude",
                  "Concentration of the total on one item" = "concentration")),
    conditionalPanel(
      "input.compare_var == 'n_ni'",
      sliderInput("compare_k", "Range of non-invariant items", 1, 7,
                  c(2, 7), step = 1, ticks = FALSE)),
    conditionalPanel(
      "input.compare_var == 'n_ni_groups'",
      sliderInput("compare_g", "Range of non-invariant groups", 1, 1,
                  c(1, 1), step = 1, ticks = FALSE),
      helpText("Requires at least 3 groups (the first group always stays ",
               "invariant).")),
    conditionalPanel(
      "input.compare_var == 'p'",
      textInput("compare_p", "p values (comma-separated)", "4, 6, 8, 10, 12, 14")),
    conditionalPanel(
      "input.compare_var == 'concentration'",
      textInput("compare_conc",
                "Concentration values from 0 to 1 (comma-separated)",
                "0, 0.25, 0.5, 0.75, 1"),
      helpText("Curves are labelled by the share of the total the largest ",
               "item ends up carrying, so 0 appears as an equal share.")),
    conditionalPanel(
      "input.compare_var == 'magnitude'",
      textInput("compare_mag",
                "Magnitude values (comma-separated; these replace the magnitude above)",
                "0.05, 0.1, 0.2, 0.3"),
      helpText("Each value is the difference added in each non-invariant ",
               "group for the selected non-invariant parameter type, ",
               "interpreted per item or as a total according to the ",
               "'Magnitude specified as' setting above.")),
    conditionalPanel(
      "input.compare_var == 'n_ni' && input.mag_mode == 'total'",
      helpText("With a total magnitude, varying the number of ",
               "non-invariant items holds the total difference constant ",
               "and splits it across more items.")),
    tags$label("Tests to plot", class = "control-label"),
    checkboxInput("include_metric", "Metric invariance (loadings)", TRUE),
    checkboxInput("include_scalar", "Scalar invariance (intercepts)", TRUE),
    checkboxInput("include_strict", "Strict invariance (residual variances)",
                  FALSE),
    sliderInput("n_range", "N per group (power computed in increments of 10)",
                50, 2000, c(50, 1000), step = 50, ticks = FALSE),
    actionButton("go", "Compute power", class = "btn-primary w-100")
  ),

  navset_card_tab(
    # panels size to their content instead of sharing the card height, which
    # would squash the tables and let the next heading overlap them
    wrapper = function(...) card_body(..., fillable = FALSE),
    nav_panel(
      "Power curves",
      uiOutput("power_plot_ui"),
      div(class = "mt-2 mb-3",
          downloadButton("dl_csv", "Results (CSV)",
                         class = "btn-sm btn-outline-secondary"),
          downloadButton("dl_png", "Plot (PNG)",
                         class = "btn-sm btn-outline-secondary"))
    ),
    nav_panel(
      "Details",
      div(class = "mb-4", h5("Scale-level tests"), DTOutput("details_scale")),
      div(class = "mb-4", h5("Item-level tests"), DTOutput("details_item")),
      helpText(class = "mt-3",
               "RMSEA (configural): population RMSEA of the unconstrained ",
               "configural model (0 whenever each group's population follows ",
               "a one-factor model). RMSEA (model): population RMSEA of the ",
               "constrained model (lavaan multi-group convention). RMSEA_d: ",
               "nested-difference RMSEA comparing the constrained model with ",
               "the model one invariance level lower, ",
               "sqrt(G x (F0_model - F0_baseline) / (df_model - df_baseline)), ",
               "scaled by the number of groups G so it is on the same metric ",
               "as the model RMSEAs. df (test): degrees of freedom of the ",
               "chi-square difference test. N for 60% / 80% power: smallest ",
               "per-group N on the computed grid reaching that power; ",
               "'> max' means it is not reached within the selected N range.")
    ),
    nav_panel(
      "Population",
      div(class = "mb-4",
          h5("Conditions in the last run"),
          helpText("Between-group difference in each non-invariant group ",
                   "(that group minus group 1), item by item, for every curve ",
                   "that was computed. Only parameters that differ somewhere ",
                   "are listed."),
          DTOutput("tbl_conditions")),
      tags$hr(),
      div(class = "mb-4",
          h5("Population values for the current sidebar settings"),
          helpText("Parameter values per item (rows) and group (columns). ",
                   "These follow the sidebar and show one setting only, not ",
                   "the compared curves.")),
      div(class = "mb-4", h5("Factor loadings"),
          div(class = "narrow-table", DTOutput("tbl_loadings"))),
      div(class = "mb-4", h5("Intercepts"),
          div(class = "narrow-table", DTOutput("tbl_intercepts"))),
      div(class = "mb-4", h5("Residual variances"),
          div(class = "narrow-table", DTOutput("tbl_residuals")))
    ),
    nav_panel(
      "About",
      markdown("
**What this app computes.** Analytic (asymptotic) power of chi-square
difference tests of measurement invariance in multiple-group one-factor
CFA, following MacCallum, Browne & Cai (2006): analysis models are fit to
the population moments, the population minimum discrepancy yields the
noncentrality parameter `nonc = (N_total - G) * (F0_model - F0_baseline)`,
and power is the tail area of the noncentral chi-square distribution beyond
the critical value of the central chi-square distribution. All tests use
significance level alpha = .05.

**Scale-level vs item-level.** *Item-level* tests constrain every loading /
intercept / residual variance to equality across groups. *Scale-level*
tests constrain only the parameter **totals** (one degree of freedom per
additional group), which is sufficient for the comparability of
unit-weighted composite scores across groups. Between-group differences
that cancel in the parameter total (the 'totals cancel out' direction
settings) leave composite-score comparisons unaffected, and tests of
scale-level invariance are insensitive to them by design.

The two approaches also respond to different features of the same total
difference. Holding the total fixed and changing how it is spread over the
items (the concentration setting) leaves scale-level power essentially
unchanged, because only the total enters its constraint, while item-level
power rises steadily as the total is concentrated on fewer items.

Each test is compared against the model one invariance level lower
(configural, then metric, then scalar, then strict), fit with the same
approach. In the plots, the solid line is the scale-level test and the
dashed line is the item-level test; when several values are compared,
color distinguishes the values.
")
    )
  )
)

server <- function(input, output, session) {

  ## keep item-count-dependent sliders in range --------------------------------
  # Raise/lower the maxima with p and G, but only ever send a value that is a
  # valid slider value: a malformed one can corrupt a range slider.
  clamp_slider <- function(id, cur, mx) {
    if (length(cur) == 2)
      updateSliderInput(session, id, max = mx, value = pmin(cur, mx))
    else if (length(cur) == 1)
      updateSliderInput(session, id, max = mx, value = min(cur, mx))
    else
      updateSliderInput(session, id, max = mx)
  }
  observeEvent(input$p, {
    mx <- max(1, input$p - 1)
    clamp_slider("n_ni", input$n_ni, mx)
    clamp_slider("compare_k", input$compare_k, mx)
  })
  observeEvent(input$n_groups, {
    mx <- max(1, input$n_groups - 1)
    clamp_slider("n_ni_groups", input$n_ni_groups, mx)
    clamp_slider("compare_g", input$compare_g, mx)
  })
  # make sure the test matching the manipulated parameters is switched on
  # which parameter classes the sidebar has ticked, in the canonical order
  ni_classes <- reactive(
    MI_PARAM_CLASSES[c(isTRUE(input$ni_loadings), isTRUE(input$ni_intercepts),
                       isTRUE(input$ni_residuals))])
  # the magnitude belonging to each class, NA when its box has been cleared
  ni_magnitudes <- reactive({
    one <- function(id) {
      v <- input[[id]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    }
    c(loadings = one("magnitude"), intercepts = one("magnitude_int"),
      residuals = one("magnitude_res"))
  })
  # ticking a parameter class switches on the test that can detect it
  observeEvent(input$ni_loadings, if (isTRUE(input$ni_loadings))
    updateCheckboxInput(session, "include_metric", value = TRUE))
  observeEvent(input$ni_intercepts, if (isTRUE(input$ni_intercepts))
    updateCheckboxInput(session, "include_scalar", value = TRUE))
  observeEvent(input$ni_residuals, if (isTRUE(input$ni_residuals))
    updateCheckboxInput(session, "include_strict", value = TRUE))

  structured_args <- reactive(list(
    p = input$p, n_groups = input$n_groups,
    base_loading = input$base_loading, base_residual = input$base_residual,
    type = ni_classes(), n_ni = input$n_ni,
    magnitude = ni_magnitudes()[["loadings"]],
    magnitude_intercept = ni_magnitudes()[["intercepts"]],
    magnitude_residual = ni_magnitudes()[["residuals"]],
    magnitude_is_total = input$mag_mode == "total",
    pattern = input$pattern,
    n_ni_groups = min(input$n_ni_groups, input$n_groups - 1),
    concentration = if (input$pattern == "uniform" && !is.null(input$conc))
      input$conc else 0))

  ## population tables ---------------------------------------------------------
  # the single population the sidebar currently describes; settings that are
  # impossible (e.g. a negative residual variance) fall back to the invariant
  # base so the tables still show something while the sidebar is mid-edit
  population <- reactive({
    tryCatch(do.call(build_structured_params, structured_args()),
             error = function(e) default_params(input$p, input$n_groups,
                                                input$base_loading))
  })

  param_table <- function(what) {
    renderDT({
      datatable(round(population()[[what]], 4), rownames = TRUE,
                options = list(dom = "t", paging = FALSE, ordering = FALSE,
                               autoWidth = FALSE,
                               columnDefs = list(list(className = "dt-right",
                                                      targets = "_all"))))
    })
  }
  output$tbl_loadings   <- param_table("loadings")
  output$tbl_intercepts <- param_table("intercepts")
  output$tbl_residuals  <- param_table("residuals")

  ## main computation ----------------------------------------------------------
  # Build one population per compared value, skipping values that are
  # impossible for the current settings (e.g. a single non-invariant item
  # with a cancelling direction) instead of losing the whole comparison.
  build_conditions <- function(args, var, values) {
    if (var == "none" || !length(values))
      return(list(conditions = structured_sweep(args, "none", NULL),
                  skipped = character(0)))
    out <- list(); skipped <- character(0)
    for (v in values) {
      one <- tryCatch(structured_sweep(args, var, v), error = function(e) NULL)
      if (is.null(one)) skipped <- c(skipped, format(v)) else out[names(one)] <- one
    }
    list(conditions = out, skipped = skipped)
  }

  res_rv <- reactiveVal(NULL)
  results <- reactive(res_rv())

  observeEvent(input$go, {
    req(input$base_loading, input$base_residual)
    res_rv(NULL)                       # never leave a stale plot on screen
    tests <- c(if (isTRUE(input$include_metric)) "metric",
               if (isTRUE(input$include_scalar)) "scalar",
               if (isTRUE(input$include_strict)) "strict")
    if (!length(tests)) {
      showNotification("Select at least one test to plot.", type = "error",
                       duration = 12)
      return()
    }
    # only the ticked classes need a magnitude; a magnitude sweep supplies
    # its own values, so the sidebar boxes are hidden and irrelevant then
    if (input$compare_var != "magnitude" &&
        anyNA(ni_magnitudes()[ni_classes()])) {
      showNotification(paste("Enter a magnitude for every ticked non-invariant",
                             "parameter type."), type = "error", duration = 12)
      return()
    }
    n_seq <- seq(input$n_range[1], input$n_range[2], by = 10)

    built <- tryCatch({
      values <- switch(input$compare_var,
                       none = NULL,
                       n_ni = seq(min(input$compare_k), max(input$compare_k)),
                       n_ni_groups = seq(min(input$compare_g), max(input$compare_g)),
                       p    = num_list(input$compare_p),
                       magnitude = num_list(input$compare_mag),
                       concentration = num_list(input$compare_conc))
      if (input$compare_var == "concentration")
        values <- values[values >= 0 & values <= 1]
      if (input$compare_var == "p")
        values <- values[values >= 3 & values > input$n_ni]
      if (input$compare_var == "n_ni_groups")
        values <- values[values >= 1 & values <= input$n_groups - 1]
      if (input$compare_var != "none" && !length(values))
        stop("No usable values to compare. Check the range or list of values.")
      if (input$compare_var != "none" && length(values) == 1)
        showNotification(
          if (input$compare_var == "n_ni_groups")
            paste("Only one value to compare: with 2 groups just one group",
                  "can be non-invariant. Comparing the number of",
                  "non-invariant groups needs at least 3 groups.")
          else "Only one value to compare, so a single curve is drawn.",
          type = "warning", duration = 12)
      b <- build_conditions(structured_args(), input$compare_var, values)
      if (!length(b$conditions))
        stop(paste("None of the compared values work with the current",
                   "settings (skipped:", paste(b$skipped, collapse = ", "), ")."))
      b
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 12)
      NULL
    })
    req(built)

    if (length(built$skipped))
      showNotification(paste0("Skipped ", length(built$skipped),
                              " value(s) impossible with the current settings: ",
                              paste(built$skipped, collapse = ", "),
                              ". Cancelling directions need at least 2 ",
                              "non-invariant items."),
                       type = "warning", duration = 12)

    param_list <- built$conditions
    out <- withProgress(message = "Fitting population models...", value = 0, {
      tryCatch(
        mi_power_sweep(param_list, tests = tests, n_seq = n_seq,
                       alpha = .05,
                       .progress = function(i, nm)
                         incProgress(1 / length(param_list), detail = nm)),
        error = function(e) {
          showNotification(paste("Computation failed:", conditionMessage(e)),
                           type = "error", duration = 12)
          NULL
        })
    })
    if (!is.null(out)) out$params <- param_list
    res_rv(out)
  })

  n_tests <- reactive({
    req(results())
    length(unique(results()$curves$test))
  })

  output$power_plot_ui <- renderUI({
    req(results())
    plotOutput("power_plot",
               height = if (n_tests() > 2) "840px" else "470px")
  })
  legend_title <- reactive({
    ttl <- COMPARE_TITLES[input$compare_var]
    if (is.na(ttl)) NULL else unname(ttl)
  })

  output$power_plot <- renderPlot({
    req(results())
    plot_power(results()$curves, ncol = min(2, n_tests()),
               color_title = legend_title())
  }, res = 100)

  # Table contents are built by plain functions so they can be checked
  # directly; the render functions only wrap them for display.
  details_data <- function(appr) {
    req(results())
    d <- results()$details
    d <- d[d$approach == appr, , drop = FALSE]
    n_max <- max(results()$curves$n)
    out <- data.frame(
      Test = unname(TEST_NAMES[d$test]),
      `RMSEA (configural)` = round(d$rmsea_configural, 4),
      `RMSEA (model)` = round(d$rmsea, 4),
      RMSEA_d = round(d$rmsea_d, 4),
      `df (test)` = d$df_test,
      `N for 60% power` = ifelse(is.na(d$n_for_60),
                                 paste0("> ", n_max), as.character(d$n_for_60)),
      `N for 80% power` = ifelse(is.na(d$n_for_80),
                                 paste0("> ", n_max), as.character(d$n_for_80)),
      check.names = FALSE)
    # keep the label visible even when a comparison collapsed to one curve,
    # so it is clear which value was actually computed
    if (any(d$condition != "current"))
      out <- cbind(Curve = d$condition, out)
    out
  }

  details_table <- function(appr) {
    renderDT(datatable(details_data(appr), rownames = FALSE,
                       options = list(dom = "t", paging = FALSE,
                                      ordering = FALSE)))
  }
  output$details_scale <- details_table("scale")
  output$details_item  <- details_table("item")

  # item-by-item differences behind each curve, so the compared conditions
  # can be checked against what was intended
  PARAM_NAMES <- c(loadings = "Loading", intercepts = "Intercept",
                   residuals = "Residual variance")

  conditions_data <- function() {
    pl <- results()$params
    if (is.null(pl)) return(NULL)
    # curves can have different numbers of items (comparing over p), so lay
    # every row out on the union of the item names, longest curve first; items
    # a curve does not have stay empty.
    per_curve <- lapply(pl, function(prm) rownames(prm$loadings))
    items <- unique(unlist(per_curve[order(-lengths(per_curve))],
                           use.names = FALSE))
    rows <- list()
    for (type in names(PARAM_NAMES)) {
      diffs <- lapply(pl, function(prm)
        round(prm[[type]][, prm$n_groups] - prm[[type]][, 1], 4))
      if (all(vapply(diffs, function(d) all(d == 0), logical(1)))) next
      for (i in seq_along(diffs)) {
        d <- setNames(rep(NA_real_, length(items)), items)
        d[names(diffs[[i]])] <- diffs[[i]]
        rows[[length(rows) + 1]] <- data.frame(
          Parameter = PARAM_NAMES[[type]], Curve = names(pl)[i],
          as.list(d),
          Total = round(sum(diffs[[i]]), 4), check.names = FALSE)
      }
    }
    if (!length(rows)) return(NULL)
    out <- do.call(rbind, rows)
    # a lone population is named "current", which labels nothing: drop the
    # column, exactly as details_data() does
    if (all(out$Curve == "current")) out$Curve <- NULL
    out
  }

  output$tbl_conditions <- renderDT({
    validate(need(!is.null(results()$params),
                  "Press 'Compute power' to see the conditions."))
    out <- conditions_data()
    validate(need(!is.null(out),
                  "All groups have identical parameters in every curve."))
    datatable(out, rownames = FALSE,
              options = list(dom = "t", paging = FALSE, ordering = FALSE,
                             scrollX = TRUE,
                             columnDefs = list(list(className = "dt-right",
                                                    targets = "_all"))))
  })

  output$dl_csv <- downloadHandler(
    filename = function() "mi_power_curves.csv",
    content = function(file) write.csv(results()$curves, file, row.names = FALSE))

  output$dl_png <- downloadHandler(
    filename = function() "mi_power_plot.png",
    content = function(file) {
      ggsave(file, plot_power(results()$curves, ncol = min(2, n_tests()),
                              color_title = legend_title()),
             width = if (n_tests() > 1) 11 else 7,
             height = if (n_tests() > 2) 9 else 5, dpi = 300)
    })
}

shinyApp(ui, server)
