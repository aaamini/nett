# function to compute the rejection rates
get_reject_freq <- function(Tstat, th_vec, twosided=T) {
  # Tstat is a vector, th could be a vector of thresholds threshold
  if (twosided) Tstat <- abs(Tstat)
  sapply(th_vec, function(th) mean(Tstat > th))
}

# function to compute the ROC
get_roc <- function(T0, T1, twosided=T) {
  T0_sorted <- sort(unique(T0), decreasing = T)
  data.frame(FPR = get_reject_freq(T0, T0_sorted, twosided = twosided),
             TPR = get_reject_freq(T1, T0_sorted, twosided = twosided))
}

# methods is a named list of functions
# gen_null_data and gen_null_data are functions that generate data under the null and alternative
#' @export
simulate_roc = function(methods, gen_null_data, gen_alt_data, nruns = 100,
                        twosided = T,
                        core_count = parallel::detectCores() - 1) {
  # method_names <- gsub("p","+",names(methods))
  method_names <- names(methods)
  n_methods <- length(methods)
  if (length(twosided) == 1) {
    twosided = rep(twosided, n_methods)
  }

  simulate_roc_run = function(j) {
    null_data = gen_null_data()
    alt_data = gen_alt_data()

    result = NULL
    for (r in 1:n_methods) {
      mtd_name = method_names[r]
      mtd_fun = methods[[r]]
      result = dplyr::bind_rows(result, tibble::tibble(method = mtd_name, H=0, tstat=mtd_fun(null_data)) )
      result = dplyr::bind_rows(result, tibble::tibble(method = mtd_name, H=1, tstat=mtd_fun(alt_data)) )
    }
    result
  }

  elapsed_time <- as.numeric(
    system.time(
      result <- do.call(rbind, parallel::mclapply(1:nruns, simulate_roc_run, mc.cores = core_count))
    )["elapsed"]
  )

  roc_results <- NULL
  for (r in 1:n_methods) {
    mtd_name <- method_names[r]
    curr_result <- dplyr::filter(result, method == mtd_name)
    T0 <- dplyr::filter(curr_result, H==0)$tstat
    T1 <- dplyr::filter(curr_result, H==1)$tstat

    roc_results <- rbind(roc_results,
                         get_roc(T0, T1, twosided = twosided[r]) %>%
                           tibble::add_column(method=mtd_name) )
  }
  list(roc=roc_results, raw=result, elapsed_time=elapsed_time)
}

#' @export
plot_roc <- function(roc_results, method_names=NULL) {
  if (!is.null(method_names)){
    roc_results = roc_results %>%
      dplyr::mutate(method = factor(method, levels = method_names))
  } else {
    roc_results = roc_results %>%
      dplyr::mutate(method = factor(method))
  }
  p = roc_results %>%
    ggplot2::ggplot(ggplot2::aes(x = FPR, y = TPR, color = method)) +
    ggplot2::geom_line(size=1.1)  +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::coord_fixed() +
    ggplot2::geom_abline(intercept =0 , slope = 1, linetype="dashed") +
    ggplot2::scale_x_continuous(limits = c(0,1.01), expand = c(0,0)) +
    ggplot2::scale_y_continuous(limits = c(0,1.01),expand = c(0,0)) +
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = c(0.85, 0.2),
      legend.text = ggplot2::element_text(size=12),
      text = ggplot2::element_text(size=16)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
      keywidth = 0.25,
      keyheight = 0.25,
      default.unit = "inch")
    )
  p
}
