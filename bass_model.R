# Bass fitting on annual cumulative counts.
# The fitting formula, starts, bounds and retry grid follow the study analysis.

fit_bass <- function(series) {
  cumulative <- series$cumulative_sum
  previous <- c(0, head(cumulative, -1L))
  ols_data <- data.frame(A_t = cumulative - previous,
                         N_t_minus_1 = previous)
  ols_data <- ols_data[series$time_model >= 1, , drop = FALSE]
  ols_fit <- lm(A_t ~ N_t_minus_1 + I(N_t_minus_1^2), data = ols_data)
  beta <- coef(ols_fit)
  discriminant <- beta[2]^2 - 4 * beta[1] * beta[3]

  if (discriminant > 0) {
    roots <- c((-beta[2] + sqrt(discriminant)) / (2 * beta[3]),
               (-beta[2] - sqrt(discriminant)) / (2 * beta[3]))
    candidates <- roots[roots > 0]
    m_start <- if (length(candidates)) min(candidates) else 1.5 * max(cumulative)
    p_start <- beta[1] / m_start
    q_start <- beta[2] + p_start
    p_start <- max(p_start, 0.0001)
    q_start <- max(q_start, 0.0001)
  } else {
    m_start <- 1.2 * max(cumulative)
    p_start <- 0.001
    q_start <- 0.3
  }

  # t = 1 is the first annual cumulative observation; N(0) = 0.
  bass_formula <- N_t ~ m * (1 - exp(-(p + q) * tempo)) /
    (1 + (q / p) * exp(-(p + q) * tempo))
  lower <- c(m = max(cumulative), p = 0, q = 0)
  fit_once <- function(start) {
    tryCatch(minpack.lm::nlsLM(
      formula = bass_formula,
      data = list(N_t = cumulative, tempo = series$time_model),
      start = start, lower = lower,
      control = minpack.lm::nls.lm.control(maxiter = 1000)
    ), error = identity)
  }

  result <- fit_once(list(m = m_start, p = p_start, q = q_start))
  if (!inherits(result, "error")) return(result)

  grid <- expand.grid(m_mult = c(0.7, 0.9, 1, 1.1, 1.3),
                      p_mult = c(0.5, 1, 1.5, 2),
                      q_mult = c(0.5, 1, 1.5, 2))
  for (i in seq_len(nrow(grid))) {
    result <- fit_once(list(
      m = max(m_start * grid$m_mult[i], 1.5 * max(cumulative)),
      p = max(p_start * grid$p_mult[i], 0.00001),
      q = max(q_start * grid$q_mult[i], 0.00001)))
    if (!inherits(result, "error")) return(result)
  }
  stop("No possible fitting: ", conditionMessage(result))
}

bass_intervals <- function(model) {
  pars <- coef(model)
  result <- c(se_m = NA_real_, se_p = NA_real_, se_q = NA_real_,
    m_ci_low = NA_real_, m_ci_high = NA_real_,
    p_ci_low = NA_real_, p_ci_high = NA_real_, q_ci_low = NA_real_,
    q_ci_high = NA_real_)
  covariance <- tryCatch(vcov(model), error = function(e) NULL)
  if (is.null(covariance) || any(!is.finite(covariance))) return(result)
  result[c("se_m", "se_p", "se_q")] <- sqrt(pmax(diag(covariance)[c("m", "p", "q")], 0))
  z <- qnorm(0.975)
  for (parameter in c("m", "p", "q")) {
    if (pars[parameter] <= 0) next
    variance <- covariance[parameter, parameter] / pars[parameter]^2
    if (!is.finite(variance) || variance < 0) next
    result[paste0(parameter, "_ci_low")] <- exp(log(pars[parameter]) - z * sqrt(variance))
    result[paste0(parameter, "_ci_high")] <- exp(log(pars[parameter]) + z * sqrt(variance))
  }
  result
}
