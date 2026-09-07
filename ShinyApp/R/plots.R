# ---------------------------------------------------------------------------
# plots.R -- ggplot2 presentation of mi_power() / mi_power_sweep() results.
# Solid = scale-level, dashed = item-level.
# When several conditions are overlaid they all stay in one panel, so that
# the two approaches remain directly comparable: each approach keeps its own
# hue (blue / orange) and the compared values run dark-to-light within that
# hue, with the value printed at the right-hand tip of every curve. A
# qualitative palette was tried first and fails past about five conditions,
# because nothing tells the reader which colour is the larger value.
# ---------------------------------------------------------------------------

# CVD-validated palettes (Okabe-Ito based)
MI_APPROACH_COLORS <- c(scale = "#0072B2", item = "#D55E00")

# Dark-to-light ramps anchored on the two approach colours. The light ends
# stop well short of white: a paler line is unreadable on the panel.
MI_CONDITION_RAMPS <- list(scale = c("#052F5F", "#1F6FB4", "#8CBFE4"),
                           item  = c("#6B2200", "#D55E00", "#FBA968"))

MI_TEST_LABELS <- c(metric = "Metric invariance (loadings)",
                    scalar = "Scalar invariance (intercepts)",
                    strict = "Strict invariance (residual variances)")
MI_APPROACH_LABELS <- c(scale = "Scale-level", item = "Item-level")

#' `k` colours for one approach, darkest = first condition.
mi_condition_colors <- function(approach, k) {
  if (k <= 1) return(unname(MI_APPROACH_COLORS[[approach]]))
  grDevices::colorRampPalette(MI_CONDITION_RAMPS[[approach]])(k)
}

#' Shortest label that still identifies a condition ("p = 12" -> "12").
#'
#' Condition names come from structured_sweep() and always carry the compared
#' value; only that value is printed at the curve tip, because the full name
#' is too wide to sit inside the panel. Falls back to the full names if
#' trimming them would make two curves share a label.
mi_tip_labels <- function(x) {
  m <- regexpr("[0-9]+(\\.[0-9]+)?%?", x)
  short <- x
  short[m > 0] <- regmatches(x, m)
  if (anyDuplicated(short)) x else short
}

#' End points of every curve, nudged apart so converging curves stay legible.
#'
#' Labels are spread within each facet only, and kept inside the 0-1 panel:
#' the gap shrinks if there are too many curves to fit at the usual spacing.
mi_tip_positions <- function(curves, gap = .038) {
  d <- curves[curves$n == max(curves$n), , drop = FALSE]
  out <- lapply(split(d, d$test, drop = TRUE), function(dd) {
    if (!nrow(dd)) return(dd)
    dd <- dd[order(-dd$power), , drop = FALSE]
    g <- min(gap, .96 / max(1, nrow(dd) - 1))
    y <- dd$power
    for (i in seq_along(y)[-1])
      if (y[i] > y[i - 1] - g) y[i] <- y[i - 1] - g
    # Nudging only ever moves a label down, and a run of near-coincident
    # curves cascades, so the column can end up taller than the panel or
    # hanging below it. Squeeze it to fit, then slide it back inside.
    span <- diff(range(y))
    if (span > 1) y <- (y - min(y)) / span
    y <- y - min(0, min(y))
    y <- y - max(0, max(y) - 1)
    dd$tip_y <- pmin(pmax(y, 0), 1)             # guard the scale limits
    dd
  })
  do.call(rbind, out)
}

#' Power curves, faceted by invariance test.
#'
#' @param curves data frame from mi_power()$curves (optionally with a
#'   `condition` column from mi_power_sweep())
#' @param target horizontal reference line (default .8); NULL to omit
#' @param color_title what the compared curves vary in (e.g. "Magnitude"),
#'   used to say what the tip labels are; NULL falls back to a generic phrase
#' @param base_size base font size of the plot
plot_power <- function(curves, target = .8, ncol = NULL, color_title = NULL,
                       base_size = 11) {
  curves$test <- droplevels(factor(curves$test,
                                   levels = names(MI_TEST_LABELS),
                                   labels = MI_TEST_LABELS))
  curves$approach <- factor(curves$approach, levels = names(MI_APPROACH_LABELS),
                            labels = MI_APPROACH_LABELS)
  multi <- "condition" %in% names(curves) &&
    length(unique(curves$condition)) > 1

  if (multi) {
    curves$condition <- factor(curves$condition, levels = unique(curves$condition))
    # one series per approach x condition, ordered approach-first so the two
    # colour ramps can simply be concatenated below
    appr_lv <- levels(curves$approach)
    cond_lv <- levels(curves$condition)
    series_lv <- unlist(lapply(appr_lv, function(a) paste(a, cond_lv, sep = " | ")))
    curves$series <- factor(paste(curves$approach, curves$condition, sep = " | "),
                            levels = series_lv)
    aes_map <- ggplot2::aes(x = n, y = power, color = series,
                            linetype = approach, group = series)
  } else {
    aes_map <- ggplot2::aes(x = n, y = power, color = approach,
                            linetype = approach, group = approach)
  }

  # start the axis exactly at the smallest N and always label it, instead of
  # leaving unlabelled padding before the first curve
  n_rng <- range(curves$n)
  x_breaks <- pretty(n_rng, n = 5)
  x_breaks <- x_breaks[x_breaks > n_rng[1] & x_breaks <= n_rng[2]]
  if (diff(n_rng) > 0)                       # drop a break that would crowd it
    x_breaks <- x_breaks[(x_breaks - n_rng[1]) / diff(n_rng) > .08]
  x_breaks <- sort(unique(c(n_rng[1], x_breaks)))

  g <- ggplot2::ggplot(curves, aes_map)
  if (!is.null(target))
    g <- g + ggplot2::geom_hline(yintercept = target, linewidth = .4,
                                 linetype = "dotted", color = "grey55")
  g <- g + ggplot2::geom_line(linewidth = .9)
  if (multi) {
    tips <- mi_tip_positions(curves)
    # trim on the distinct conditions: every condition appears once per
    # approach, and duplicates there must not defeat the ambiguity check
    tip_of <- setNames(mi_tip_labels(cond_lv), cond_lv)
    tips$tip_label <- unname(tip_of[as.character(tips$condition)])
    g <- g + ggplot2::geom_text(
      data = tips,
      mapping = ggplot2::aes(x = n, y = tip_y, label = tip_label, color = series),
      hjust = -.3, size = base_size / 3.4, fontface = "bold",
      show.legend = FALSE, inherit.aes = FALSE)
  }
  g <- g +
    ggplot2::facet_wrap(~test, ncol = ncol) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, .2), expand = c(0.01, 0)) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks, limits = n_rng,
      # compared curves are labelled at their tips and need room for it
      expand = ggplot2::expansion(mult = c(0, if (multi) .10 else .03))) +
    ggplot2::scale_linetype_manual(values = c("Scale-level" = "solid",
                                              "Item-level" = "42"),
                                   name = if (multi) "Approach" else NULL) +
    ggplot2::labs(x = "N per group", y = "Power") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      strip.text = ggplot2::element_text(face = "bold"),
      # the axis starts flush at the smallest N, so neighbouring panels need
      # room for their first and last labels
      panel.spacing.x = grid::unit(1.6, "lines"),
      panel.spacing.y = grid::unit(1, "lines"),
      # only the approach key is left in the legend; the compared values are
      # named at the curve tips instead
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "top",
      legend.key.width = grid::unit(2.2, "lines"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.spacing.y = grid::unit(.4, "lines"))

  if (multi) {
    key_of <- setNames(names(MI_APPROACH_LABELS), unname(MI_APPROACH_LABELS))
    cols <- unlist(lapply(appr_lv,
                          function(a) mi_condition_colors(key_of[[a]], length(cond_lv))))
    names(cols) <- series_lv
    tip <- mi_tip_labels(cond_lv)
    # the colour key would need one row per approach x condition; the tip
    # labels already name every curve, so only the ramp itself is explained
    g <- g + ggplot2::scale_color_manual(values = cols, name = color_title,
                                         guide = "none") +
      ggplot2::guides(linetype = ggplot2::guide_legend(order = 1)) +
      ggplot2::labs(caption = sprintf(
        "Curve tips give %s. Blue = scale-level, orange = item-level; within each, dark = %s to light = %s.",
        if (is.null(color_title)) "the compared value" else tolower(color_title),
        tip[1], tip[length(tip)])) +
      ggplot2::theme(plot.caption = ggplot2::element_text(
        hjust = 0, color = "grey35", size = base_size * .8,
        margin = ggplot2::margin(t = 8)))
  } else {
    g <- g + ggplot2::scale_color_manual(
      values = setNames(unname(MI_APPROACH_COLORS[c("scale", "item")]),
                        unname(MI_APPROACH_LABELS[c("scale", "item")])),
      name = NULL)
  }
  g
}
