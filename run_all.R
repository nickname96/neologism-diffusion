# Run from a terminal: Rscript run_all.R
# Or, in R/RStudio, set the working directory to this folder and source this file.
# Packages: minpack.lm, ggplot2, ggrepel, patchwork, paletteer, scales.

packages <- c("minpack.lm", "ggplot2", "ggrepel", "patchwork", "paletteer", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the missing packages first: ", paste(missing, collapse = ", "))
}
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)))
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
release_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/"))
} else getwd()
if (!file.exists(file.path(release_dir, "bass_model.R"))) {
  stop("In R/RStudio, set the working directory to the release folder first.")
}
source(file.path(release_dir, "bass_model.R"))
library(ggplot2)
library(patchwork)

# Default settings used for the study's count-based Bass analysis.
settings <- list(first_year = 2006L, last_year = 2024L,
                 min_total_uses = 500L, min_observation_years = 6L)
result_dir <- file.path(release_dir, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, filename) {
  write.csv(x, file.path(result_dir, filename), row.names = FALSE,
            na = "", fileEncoding = "UTF-8")
}
read_counts <- function(filename, platform) {
  x <- read.csv(file.path(release_dir, "data", filename),
                stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  stopifnot(all(c("year", "neologism", "count") %in% names(x)),
            all(is.finite(x$year)), all(is.finite(x$count)),
            all(x$count >= 0), all(x$count == floor(x$count)),
            !anyNA(x$neologism), all(nzchar(x$neologism)),
            !anyDuplicated(x[c("neologism", "year")]))
  x <- x[x$year >= settings$first_year & x$year <= settings$last_year, ]
  x$platform <- platform
  x
}
annual <- rbind(read_counts("time_series_neologisms_reddit.csv", "Reddit"),
                read_counts("time_series_neologisms_youtube.csv", "YouTube"))
categories <- read.csv(file.path(release_dir, "data", "neologism_categories.csv"),
                       stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
stopifnot(!anyDuplicated(categories$neologism),
          all(annual$neologism %in% categories$neologism))

trajectories <- parameters <- curves <- list()
numeric_fields <- c("m", "p", "q", "se_m", "se_p", "se_q",
  "m_ci_low", "m_ci_high", "p_ci_low", "p_ci_high", "q_ci_low", "q_ci_high",
  "rmse", "nrmse", "r_squared")

for (platform in c("Reddit", "YouTube")) {
  words <- sort(unique(annual$neologism[annual$platform == platform]))
  for (word in words) {
    x <- annual[annual$platform == platform & annual$neologism == word, ]
    if (!any(x$count > 0)) next
    first <- min(x$year[x$count > 0])
    years <- seq.int(first, settings$last_year)
    counts <- x$count[match(years, x$year)]
    counts[is.na(counts)] <- 0
    series <- data.frame(year = years, count = counts,
                          time_model = seq_along(years), cumulative_sum = cumsum(counts))
    total <- sum(counts)
    eligible <- total >= settings$min_total_uses &&
      length(years) >= settings$min_observation_years
    metadata <- data.frame(platform = platform, neologism = word,
      first_appearance = first,
      n_obs = length(years), total_uses = total, eligible = eligible,
      status = if (eligible) "pending" else "below_threshold", error = "")
    id <- paste(platform, word, sep = "|")

    if (eligible) {
      fitted <- tryCatch(fit_bass(series), error = identity)
      metadata$status <- if (inherits(fitted, "error")) "error" else "fitted"
      if (inherits(fitted, "error")) metadata$error <- conditionMessage(fitted)
      row <- metadata
      for (name in numeric_fields) row[[name]] <- NA_real_
      if (!inherits(fitted, "error")) {
        pars <- coef(fitted)
        row[c("m", "p", "q")] <- as.list(pars[c("m", "p", "q")])
        intervals <- bass_intervals(fitted)
        row[names(intervals)] <- as.list(intervals)
        row$rmse <- sqrt(mean(residuals(fitted)^2))
        row$nrmse <- row$rmse / diff(range(series$cumulative_sum))
        tss <- sum((series$cumulative_sum - mean(series$cumulative_sum))^2)
        row$r_squared <- if (tss > 0) 1 - sum(residuals(fitted)^2) / tss else NA_real_
        series$predicted <- as.numeric(predict(fitted, newdata = data.frame(tempo = series$time_model)))
        series$platform <- platform
        series$neologism <- word
        curves[[id]] <- series[c("platform", "neologism", "year", "time_model", "count", "cumulative_sum", "predicted")]
      }
      parameters[[id]] <- row
    }
    trajectories[[id]] <- metadata
  }
  message("Finished Bass fitting for ", platform)
}
trajectory_summary <- do.call(rbind, trajectories)
if (!length(parameters)) stop("No trajectories meet the configured inclusion requirements.")
bass_parameters <- do.call(rbind, parameters)
if (!length(curves)) stop("No model returned fitted values. Check the error column in the parameter output.")
fitted_curves <- do.call(rbind, curves)
write_result(trajectory_summary, "trajectory_summary.csv")
write_result(bass_parameters, "bass_parameters.csv")

fitted_parameters <- bass_parameters[bass_parameters$status == "fitted", ]
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
platform_summary <- do.call(rbind, lapply(c("Reddit", "YouTube"), function(platform) {
  all <- trajectory_summary[trajectory_summary$platform == platform, ]
  p <- fitted_parameters[fitted_parameters$platform == platform, ]
  data.frame(platform = platform, input_terms = nrow(all),
    eligible_trajectories = sum(all$eligible), fitted_trajectories = nrow(p),
    fitting_errors = sum(all$status == "error"),
    median_p = safe_median(p$p), median_q = safe_median(p$q),
    median_NRMSE = safe_median(p$nrmse))
}))
write_result(platform_summary, "platform_summary.csv")

# Paper-style figures.
examples <- data.frame(
  platform = c("Reddit", "YouTube", "YouTube", "YouTube"),
  neologism = c("Brain rot", "Baizuo", "McSpicy", "Sigma male"),
  caption = c("Low innovation, high imitation", "High innovation, low imitation",
              "High innovation, High imitation", "Mid innovation, high imitation"))
example_plots <- lapply(seq_len(nrow(examples)), function(i) {
  ref <- examples[i, ]
  d <- fitted_curves[fitted_curves$platform == ref$platform & fitted_curves$neologism == ref$neologism, ]
  ggplot(d, aes(time_model, cumulative_sum)) +
    geom_point(size = 3.5, colour = "darkorange", alpha = 0.7) +
    geom_line(aes(y = predicted), colour = "darkblue", linewidth = 1) +
    theme_minimal() +
    theme(aspect.ratio = 0.5, text = element_text(size = 10)) +
    labs(x = NULL, y = NULL, title = paste(ref$neologism, "-", ref$platform), caption = ref$caption)
})
ggsave(file.path(figure_dir, "example_fits.png"),
       wrap_plots(example_plots, ncol = 2), width = 6, height = 4, dpi = 200, bg = "white")

plot_data <- merge(fitted_parameters, categories, by = "neologism", sort = FALSE)
plot_data$platform <- factor(plot_data$platform, levels = c("Reddit", "YouTube"))
plane_data <- plot_data[is.finite(plot_data$r_squared) & plot_data$r_squared >= 0.95 &
                         plot_data$p > 0 & plot_data$q > 0, ]
p_plane <- ggplot(plane_data, aes(p, q, colour = category, label = label_to_show)) +
  geom_point(size = 2, alpha = 0.5) +
  facet_wrap(~platform) +
  labs(x = "Innovation", y = "Imitation", colour = NULL) +
  ggrepel::geom_text_repel(show.legend = FALSE, seed = 42) +
  theme_bw() +
  theme(legend.position = "bottom", aspect.ratio = 0.5,
        strip.background = element_blank(), strip.text = element_text(face = "bold")) +
  scale_x_log10(labels = scales::label_log()) + scale_y_log10() +
  paletteer::scale_colour_paletteer_d("nationalparkcolors::Badlands") +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1)))

reddit <- fitted_parameters[fitted_parameters$platform == "Reddit", c("neologism", "p", "q")]
youtube <- fitted_parameters[fitted_parameters$platform == "YouTube", c("neologism", "p", "q")]
paired <- merge(reddit, youtube, by = "neologism", suffixes = c("_Reddit", "_YouTube"))
paired <- merge(paired, categories, by = "neologism")
parameter_comparison <- function(parameter, title) {
  d <- paired
  d$Reddit <- d[[paste0(parameter, "_Reddit")]]
  d$YouTube <- d[[paste0(parameter, "_YouTube")]]
  d <- d[is.finite(d$Reddit) & d$Reddit > 0 & is.finite(d$YouTube) & d$YouTube > 0, ]
  d$type_label <- title
  p <- ggplot(d, aes(Reddit, YouTube, colour = category)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_point(size = 2, alpha = 0.5) +
    ggrepel::geom_text_repel(aes(label = label_to_show), show.legend = FALSE, seed = 42) +
    theme_bw() + labs(x = "Reddit", y = "YouTube", colour = NULL) +
    paletteer::scale_colour_paletteer_d("nationalparkcolors::Badlands") +
    theme(aspect.ratio = 1, strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    facet_wrap(~type_label) + guides(colour = "none")
  if (parameter == "p") {
    p + scale_x_log10(labels = scales::label_log()) + scale_y_log10(labels = scales::label_log())
  } else p + scale_x_log10() + scale_y_log10()
}
parameter_figure <- (p_plane / (parameter_comparison("q", "Imitation") +
                                 parameter_comparison("p", "Innovation"))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.tag = element_text(face = "bold"), plot.tag.position = c(0.02, 0.98))
ggsave(file.path(figure_dir, "bass_parameters.png"), parameter_figure,
       width = 8, height = 8, dpi = 200, bg = "white")

error_figure <- ggplot(fitted_parameters, aes(platform, nrmse)) +
  geom_boxplot(width = 0.5, fill = "white", colour = "grey20") +
  theme_bw() + theme(aspect.ratio = 1) + labs(x = NULL, y = "NRMSE")
ggsave(file.path(figure_dir, "fit_error.png"), error_figure,
       width = 5, height = 5, dpi = 200, bg = "white")

report <- c("BASS ANALYSIS SUMMARY", "",
  sprintf("Study period: %d--%d.", settings$first_year, settings$last_year),
  sprintf("Inclusion: at least %d term uses and %d annual observations from first observed use.",
          settings$min_total_uses, settings$min_observation_years),
  "Intervening years with no recorded occurrences are included as zero counts.",
  "Cumulative counts are fitted directly; m, p and q are estimated jointly.",
  "The observed cumulative total supplies the lower bound for m.",
  "NRMSE is RMSE divided by the observed range of cumulative counts.",
  "Parameter intervals are approximate 95% log-scale covariance-based intervals.",
  "", 
  capture.output(print(platform_summary, row.names = FALSE)), "",
  "FIGURES",
  "example_fits.png: the four example trajectories used in the paper.",
  "bass_parameters.png: paper-style parameter plane and cross-platform parameter comparisons.",
  "fit_error.png: distribution of NRMSE by platform.", "",
  "Counts describe term-comment matches, not distinct users. A comment may match several terms.",
  "Fitted parameters summarize the observed curves; the figures are descriptive.", "",
  paste("R:", getRversion()),
  vapply(packages, function(p) paste(p, as.character(utils::packageVersion(p))), character(1))
)
writeLines(report, file.path(result_dir, "summary.txt"), useBytes = TRUE)
message("Done. Results: ", normalizePath(result_dir, winslash = "/"))
print(platform_summary, row.names = FALSE)
