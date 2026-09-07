# ---------------------------------------------------------------------------
# models.R -- population specification and lavaan model generation for
# scale-level vs item-level measurement invariance power analysis.
#
# Every model is generated for an arbitrary number of items (p) and
# groups (G) rather than written out as a fixed model string.
#
# A population is a `mi_params` list:
#   p, n_groups, and three p x G matrices: loadings, intercepts, residuals.
# Latent variance is fixed to 1 and latent mean to 0 in every population
# group; the observed moments are computed directly:
#   Sigma_g = lambda_g lambda_g' + diag(theta_g),  mu_g = nu_g
# ---------------------------------------------------------------------------

MI_LEVELS <- c("configural", "metric", "scalar", "strict")

#' Fully invariant population with equal parameters in all groups.
default_params <- function(p, n_groups = 2,
                           base_loading = .7,
                           base_intercept = 0,
                           base_residual = 1 - base_loading^2) {
  stopifnot(p >= 3, n_groups >= 2, base_residual > 0)
  mk <- function(x) matrix(x, nrow = p, ncol = n_groups,
                           dimnames = list(paste0("y", seq_len(p)),
                                           paste0("group", seq_len(n_groups))))
  structure(list(p = p, n_groups = n_groups,
                 loadings   = mk(base_loading),
                 intercepts = mk(base_intercept),
                 residuals  = mk(base_residual)),
            class = "mi_params")
}

#' Add non-invariance to the last `n_ni` items of one focal group.
#'
#' @param type    one of "loadings", "intercepts", "residuals"
#' @param n_ni    number of non-invariant items (item 1 is always kept
#'                invariant: it anchors the latent scale in analysis models)
#' @param delta   per-item difference added in the focal group
#' @param pattern "uniform": all items get +delta.
#'                "cancel": differences of about +/-delta spread across the
#'                items in opposite directions (alternating signs, recentered
#'                so the parameter total is exactly unchanged for any n_ni;
#'                for even n_ni this is exactly +delta, -delta, ...).
#'                "cancel_one": n_ni - 1 items get +delta and one item gets
#'                -(n_ni - 1) * delta, so the total is unchanged.
#'                "alternating" is accepted as an alias for "cancel".
#'                Cancelling patterns require n_ni >= 2.
#' @param concentration how unevenly the total difference is spread over the
#'                non-invariant items, holding that total fixed (uniform
#'                pattern only). 0 splits it evenly; 1 puts all of it on a
#'                single item; in between, the largest item's share is
#'                pi = 1/n_ni + concentration * (1 - 1/n_ni) and the rest is
#'                split evenly over the remaining items.
#' @param focal_groups which group(s) deviate (default: the last one); every
#'                focal group receives the same differences, so all focal
#'                groups are identical to each other and differ from group 1
apply_noninvariance <- function(params, type, n_ni, delta,
                                pattern = c("uniform", "cancel", "cancel_one",
                                            "alternating"),
                                focal_groups = params$n_groups,
                                concentration = 0) {
  pattern <- match.arg(pattern)
  stopifnot(length(concentration) == 1, !is.na(concentration),
            concentration >= 0, concentration <= 1)
  if (pattern == "alternating") pattern <- "cancel"
  type <- match.arg(type, c("loadings", "intercepts", "residuals"))
  stopifnot(inherits(params, "mi_params"),
            n_ni >= 0, n_ni <= params$p - 1,
            length(focal_groups) >= 1,
            all(focal_groups >= 2), all(focal_groups <= params$n_groups))
  if (n_ni == 0 || delta == 0) return(params)
  if (pattern != "uniform" && n_ni < 2)
    stop("Cancelling direction patterns require at least 2 non-invariant items.")
  items <- seq(params$p - n_ni + 1, params$p)
  deltas <- switch(pattern,
    uniform = if (concentration == 0) rep(delta, n_ni) else {
      # keep the total delta * n_ni fixed, but give one item a larger share
      total <- delta * n_ni
      share <- 1 / n_ni + concentration * (1 - 1 / n_ni)
      rest <- if (n_ni > 1) (1 - share) * total / (n_ni - 1) else numeric(0)
      c(rep(rest, n_ni - 1), share * total)
    },
    cancel = {
      s <- rep_len(c(1, -1), n_ni)
      delta * (s - mean(s))
    },
    cancel_one = delta * c(rep(1, n_ni - 1), -(n_ni - 1)))
  params[[type]][items, focal_groups] <-
    params[[type]][items, focal_groups] + deltas
  if (type == "residuals" && any(params$residuals <= 0))
    stop("Non-invariance specification produces non-positive residual variances.")
  params
}

#' Population covariance matrices and mean vectors implied by `params`.
pop_moments <- function(params) {
  stopifnot(inherits(params, "mi_params"))
  lapply(seq_len(params$n_groups), function(g) {
    lam <- params$loadings[, g]
    S <- tcrossprod(lam) + diag(params$residuals[, g])
    dimnames(S) <- list(rownames(params$loadings), rownames(params$loadings))
    list(cov = S, mean = setNames(params$intercepts[, g], rownames(S)))
  })
}

#' Plain one-factor model; item-level invariance is imposed at fit time
#' via lavaan's `group.equal`, and this model also serves as the
#' configural baseline for both approaches.
make_item_model <- function(p) {
  paste0("f1 =~ ", paste0("y", seq_len(p), collapse = " + "))
}

# label helpers: one label per item per group ---------------------------------
.lab_vec <- function(prefix, i, G) paste0(prefix, i, "_", seq_len(G))
.c_vec   <- function(x) paste0("c(", paste(x, collapse = ", "), ")")

#' Multi-group lavaan model with *scale-level* (parameter total) invariance
#' constraints up to `level` ("metric", "scalar", or "strict").
#'
#' Built for any p and G:
#'   metric:  sum of loadings equal across groups (loading 1 fixed to 1)
#'   scalar:  + first intercept anchored + intercept totals equal;
#'              latent means free in groups 2..G
#'   strict:  + residual variance totals equal
make_scale_model <- function(p, n_groups, level = c("metric", "scalar", "strict")) {
  level <- match.arg(level)
  G <- n_groups
  items <- seq_len(p)

  load_terms <- c(paste0(.c_vec(rep(1, G)), "*y1"),
                  vapply(items[-1], function(i)
                    paste0(.c_vec(.lab_vec("l", i, G)), "*y", i), character(1)))
  lines <- paste0("f1 =~ ", paste(load_terms, collapse = " + "))

  # loading-total constraints (one per non-reference group)
  lsum <- function(g) paste0("l", items[-1], "_", g, collapse = " + ")
  lines <- c(lines, vapply(2:G, function(g)
    paste(lsum(1), "==", lsum(g)), character(1)))

  if (level %in% c("scalar", "strict")) {
    lines <- c(lines,
               vapply(items, function(i)
                 paste0("y", i, " ~ ", .c_vec(.lab_vec("nu", i, G)), "*1"),
                 character(1)),
               paste0("f1 ~ ", .c_vec(c(0, rep(NA, G - 1))), "*1"))
    # anchor first intercept + intercept-total constraints
    nsum <- function(g) paste0("nu", items, "_", g, collapse = " + ")
    lines <- c(lines,
               vapply(2:G, function(g) paste0("nu1_1 == nu1_", g), character(1)),
               vapply(2:G, function(g) paste(nsum(1), "==", nsum(g)), character(1)))
  }

  if (level == "strict") {
    lines <- c(lines,
               vapply(items, function(i)
                 paste0("y", i, " ~~ ", .c_vec(.lab_vec("th", i, G)), "*y", i),
                 character(1)))
    tsum <- function(g) paste0("th", items, "_", g, collapse = " + ")
    lines <- c(lines, vapply(2:G, function(g)
      paste(tsum(1), "==", tsum(g)), character(1)))
  }

  paste(lines, collapse = "\n")
}

#' The three parameter classes that can be made non-invariant, in the order
#' they are always applied and reported.
MI_PARAM_CLASSES <- c("loadings", "intercepts", "residuals")

#' Normalise a `type` argument to a subset of MI_PARAM_CLASSES.
#'
#' Any combination of the three classes may be given, in any order and with
#' repeats; the result is de-duplicated and put in a fixed order so that two
#' equivalent requests build an identical population. NULL or character(0)
#' means "nothing is non-invariant". "both" and "all" are kept as shorthands
#' from before the classes became independently selectable.
mi_param_classes <- function(type) {
  if (is.null(type)) return(character(0))
  type <- as.character(type)
  type <- unlist(lapply(type, function(x)
    switch(x, both = c("loadings", "intercepts"), all = MI_PARAM_CLASSES, x)),
    use.names = FALSE)
  if (is.null(type)) return(character(0))
  bad <- setdiff(type, MI_PARAM_CLASSES)
  if (length(bad))
    stop("Unknown non-invariant parameter type: ", paste(bad, collapse = ", "))
  MI_PARAM_CLASSES[MI_PARAM_CLASSES %in% type]
}

#' Build a population from the structured ("preset") specification used by
#' the Shiny app: a fully invariant base, plus non-invariance in any
#' combination of the three parameter classes, on the last `n_ni` items of
#' the last `n_ni_groups` groups (each with the same differences).
#'
#' @param type which parameter classes are non-invariant: any subset of
#'   "loadings", "intercepts" and "residuals" (see mi_param_classes), or
#'   character(0) for a fully invariant population
#' @param magnitude per-item difference for the loadings, or the total
#'   across all non-invariant items if `magnitude_is_total`
#'   (Study-3-style designs)
#' @param magnitude_intercept as `magnitude`, for the intercepts; the three
#'   classes are on different scales, so each takes its own difference
#' @param magnitude_residual as `magnitude`, for the residual variances
#' @param n_ni_groups number of non-invariant groups (the last ones);
#'   at most n_groups - 1 so the reference group stays invariant
#' @param concentration 0 spreads the total difference evenly over the
#'   non-invariant items, 1 puts all of it on a single item (see
#'   apply_noninvariance); the total is the same either way. Applied to each
#'   class separately, so each class keeps its own total.
build_structured_params <- function(p, n_groups = 2,
                                    base_loading = .7,
                                    base_residual = 1 - base_loading^2,
                                    type = "loadings",
                                    n_ni = 0,
                                    magnitude = 0,
                                    magnitude_intercept = magnitude,
                                    magnitude_residual = magnitude,
                                    magnitude_is_total = FALSE,
                                    pattern = "uniform",
                                    n_ni_groups = 1,
                                    concentration = 0) {
  type <- mi_param_classes(type)
  stopifnot(n_ni_groups >= 1, n_ni_groups <= n_groups - 1)
  prm <- default_params(p, n_groups, base_loading, base_residual = base_residual)
  if (n_ni == 0 || !length(type)) return(prm)
  gs <- seq(n_groups - n_ni_groups + 1, n_groups)
  per_item <- function(m) if (magnitude_is_total) m / n_ni else m
  deltas <- c(loadings = magnitude, intercepts = magnitude_intercept,
              residuals = magnitude_residual)
  for (cls in type)
    prm <- apply_noninvariance(prm, cls, n_ni, per_item(deltas[[cls]]),
                               pattern, focal_groups = gs,
                               concentration = concentration)
  prm
}

#' Named list of populations sweeping one design factor, for overlay plots.
#'
#' @param args list of arguments for build_structured_params()
#' @param sweep_var "none", "n_ni", "p", "magnitude", "n_ni_groups", or
#'   "concentration"
structured_sweep <- function(args, sweep_var = "none", values = NULL) {
  if (sweep_var == "none" || is.null(values) || length(values) == 0) {
    return(setNames(list(do.call(build_structured_params, args)), "current"))
  }
  stopifnot(sweep_var %in% c("n_ni", "p", "magnitude", "n_ni_groups",
                             "concentration"))
  lab <- switch(sweep_var,
                n_ni = "%s non-invariant",
                p = "p = %s",
                n_ni_groups = "non-invariant groups = %s",
                concentration = "largest item = %s of total",
                magnitude = if (isTRUE(args$magnitude_is_total))
                  "total diff = %s" else "per-item diff = %s")
  # label concentrations by the share the largest item ends up carrying
  fmt <- if (sweep_var == "concentration") {
    k <- args$n_ni
    function(v) paste0(round(100 * (1 / k + v * (1 - 1 / k))), "%")
  } else format
  out <- list()
  for (v in values) {
    a <- args
    a[[sweep_var]] <- v
    # a magnitude sweep moves every parameter class that `type` includes
    if (sweep_var == "magnitude") {
      a$magnitude_intercept <- v
      a$magnitude_residual <- v
    }
    out[[sprintf(lab, fmt(v))]] <- do.call(build_structured_params, a)
  }
  out
}

#' `group.equal` argument for the item-level model at a given level.
item_group_equal <- function(level) {
  switch(level,
         configural = character(0),
         metric     = "loadings",
         scalar     = c("loadings", "intercepts"),
         strict     = c("loadings", "intercepts", "residuals"))
}
