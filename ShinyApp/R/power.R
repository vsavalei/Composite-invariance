# ---------------------------------------------------------------------------
# power.R -- analytic power of chi-square difference tests for scale-level
# and item-level measurement invariance (MacCallum, Browne & Cai, 2006).
#
# Replaces the power.plot.doer.* functions from
# "power analysis plot function.R", with three changes:
#   * computation is separated from plotting and returns tidy data frames
#   * degrees of freedom come from the fitted models, fixing the original
#     df bug in the metric-test critical value (nrow(pop.cov1-1) gave p
#     instead of p-1)
#   * the noncentrality parameter is computed from the population minimum
#     discrepancy F0 = chisq / N_big, so it is correct for any number of
#     groups regardless of lavaan's multi-group RMSEA convention:
#       nonc(N) = (N_total - G) * (F0_constrained - F0_baseline)
#     For G = 2 this is numerically identical to the manuscript's
#     n * df * rmsea^2 (verified against the published figures).
# ---------------------------------------------------------------------------

#' Fit an analysis model to the population moments with a huge sample size.
fit_pop_model <- function(model, moments, group.equal = character(0),
                          n_big = 1e7) {
  covs  <- lapply(moments, `[[`, "cov")
  means <- lapply(moments, `[[`, "mean")
  nobs  <- rep(n_big, length(moments))
  if (length(group.equal)) {
    lavaan::sem(model, sample.cov = covs, sample.mean = means,
                sample.nobs = nobs, group.equal = group.equal)
  } else {
    lavaan::sem(model, sample.cov = covs, sample.mean = means,
                sample.nobs = nobs)
  }
}

.model_stats <- function(fit, n_big, G) {
  fm <- lavaan::fitmeasures(fit, c("chisq", "df", "rmsea"))
  list(F0 = unname(fm["chisq"]) / (G * n_big),
       df = unname(fm["df"]),
       rmsea = unname(fm["rmsea"]))
}

#' Power curves for tests of measurement invariance in a given population.
#'
#' For each requested test the baseline is the model one invariance level
#' lower (configural for metric, metric for scalar, scalar for strict),
#' fit with the same approach, as in the manuscript.
#'
#' @param params  an `mi_params` population (see models.R)
#' @param tests   subset of c("metric", "scalar", "strict")
#' @param n_seq   per-group sample sizes at which to compute power
#' @param alpha   significance level of the chi-square difference test
#' @return list(curves, details): `curves` has one row per n x test x
#'   approach; `details` has population RMSEA, df and F0 per test x approach.
mi_power <- function(params,
                     tests = c("metric", "scalar"),
                     n_seq = seq(50, 1000, by = 10),
                     alpha = .05,
                     n_big = 1e7) {
  tests <- match.arg(tests, c("metric", "scalar", "strict"), several.ok = TRUE)
  tests <- MI_LEVELS[MI_LEVELS %in% tests]          # canonical order
  G <- params$n_groups
  p <- params$p
  moments <- pop_moments(params)
  item_model <- make_item_model(p)

  # fit every level needed (requested tests + their baselines), per approach
  need <- unique(unlist(lapply(tests, function(tt) {
    i <- match(tt, MI_LEVELS); MI_LEVELS[c(i - 1, i)]
  })))
  fits <- list()
  for (lev in need) {
    if (lev == "configural") {
      f <- fit_pop_model(item_model, moments, n_big = n_big)
      st <- .model_stats(f, n_big, G)
      fits[["scale.configural"]] <- fits[["item.configural"]] <- st
    } else {
      fits[[paste0("scale.", lev)]] <- .model_stats(
        fit_pop_model(make_scale_model(p, G, lev), moments, n_big = n_big),
        n_big, G)
      fits[[paste0("item.", lev)]] <- .model_stats(
        fit_pop_model(item_model, moments,
                      group.equal = item_group_equal(lev), n_big = n_big),
        n_big, G)
    }
  }

  curves <- details <- list()
  for (tt in tests) {
    base_lev <- MI_LEVELS[match(tt, MI_LEVELS) - 1]
    for (appr in c("scale", "item")) {
      con <- fits[[paste0(appr, ".", tt)]]
      bas <- fits[[paste0(appr, ".", base_lev)]]
      cfg <- fits[[paste0(appr, ".configural")]]
      df_diff <- con$df - bas$df
      dF0     <- max(con$F0 - bas$F0, 0)
      nonc    <- (G * n_seq - G) * dF0
      pow     <- pchisq(qchisq(1 - alpha, df_diff), df_diff,
                        ncp = nonc, lower.tail = FALSE)
      curves[[length(curves) + 1]] <- data.frame(
        n = n_seq, power = pow, test = tt, approach = appr)
      n_at <- function(target) {
        hit <- n_seq[which(pow >= target)[1]]
        if (is.na(hit)) NA_real_ else hit
      }
      # nested-difference RMSEA: sqrt(dF0/ddf), scaled by G so it is on the
      # same multi-group metric as lavaan's model RMSEA
      rmsea_d <- if (df_diff > 0) sqrt(G * dF0 / df_diff) else NA_real_
      details[[length(details) + 1]] <- data.frame(
        test = tt, approach = appr,
        rmsea_configural = if (is.null(cfg)) NA_real_ else cfg$rmsea,
        rmsea_baseline = bas$rmsea,
        rmsea = con$rmsea,
        rmsea_d = rmsea_d,
        df_model = con$df, df_test = df_diff,
        F0_diff = dF0, n_for_60 = n_at(.60), n_for_80 = n_at(.80))
    }
  }
  list(curves  = do.call(rbind, curves),
       details = do.call(rbind, details))
}

#' Run `mi_power` over a named list of populations (a "sweep") and combine.
#'
#' @param param_list named list of `mi_params`; names become the condition
#'   labels used for coloring in the plots.
#' @param .progress optional function(i, name) called before each condition
#'   (used by the Shiny app for its progress bar).
mi_power_sweep <- function(param_list, ..., .progress = NULL) {
  stopifnot(length(param_list) >= 1, !is.null(names(param_list)))
  out <- lapply(seq_along(param_list), function(i) {
    if (!is.null(.progress)) .progress(i, names(param_list)[i])
    res <- mi_power(param_list[[i]], ...)
    res$curves$condition  <- names(param_list)[i]
    res$details$condition <- names(param_list)[i]
    res
  })
  list(curves  = do.call(rbind, lapply(out, `[[`, "curves")),
       details = do.call(rbind, lapply(out, `[[`, "details")))
}
