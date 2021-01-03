# function to compute the rejection rates
get_reject_freq <- function(Tstat, th_vec, twosided = T) {
  # if (twosided) cat('twosided\n') else cat('onesided\n')
  # Tstat is a vector, th could be a vector of thresholds threshold
  if (twosided) Tstat <- abs(Tstat)
  sapply(th_vec, function(th) mean(Tstat > th))
}

# function to compute the ROC
get_roc <- function(T0, T1, twosided = T) {
  T0_sorted <- sort(unique(T0), decreasing = T)
  data.frame(FPR = get_reject_freq(T0, T0_sorted, twosided = twosided),
             TPR = get_reject_freq(T1, T0_sorted, twosided = twosided))
}

#' simulate data to draw ROC curves
#' @param apply_methods a function that returns a data.frame with columns "method", "tstat" and "twosided"
#' @param  gen_null_data a function that generate data under the null model
#' @param  gen_alt_data a function that generate data under the alternative model
#' @param nruns number of simulated data from the null/alternative model
#' @param core_count number of cores used in parallel computing
#' @param seed seed for random simulation
#' @return a list of result
#' \item{roc}{A data frame used to plot ROC curves with columns: method, whether a two sided test,
#' false positive rate (FPR), and true positive rate (TPR)}
#' \item{raw}{A data frame containing raw output from null and alternative models with columns:
#' method, statistics value, whether a two sided test, and the type of hypothesis}
#' \item{elapsed_time}{symstem elapsed time for generating ROC data}
#' @keywords plotting
#' @export
simulate_roc = function(apply_methods, gen_null_data, gen_alt_data,
                        nruns = 100,
                        core_count = parallel::detectCores() - 1,
                        seed = NULL) {

  simulate_roc_run = function(j) {
    null_data = gen_null_data()
    alt_data = gen_alt_data()

    dplyr::bind_rows(
      apply_methods(null_data) %>%  tibble::add_column(H=0),
      apply_methods(alt_data) %>%  tibble::add_column(H=1)
    )
  }

  if (!is.null(seed)) set.seed(seed)

  elapsed_time <- as.numeric(
    system.time(
      result <- do.call(rbind, parallel::mclapply(1:nruns, simulate_roc_run, mc.cores = core_count))
    )["elapsed"]
  )

  roc_results = result %>%
                  dplyr::group_by(method, H) %>%
                  tidyr::nest(tstat_all = tstat)  %>%
                  tidyr::pivot_wider(names_from = H, values_from = tstat_all, names_prefix="H") %>%
                  dplyr::mutate(res = purrr::map2(H0, H1, ~get_roc(unlist(.x), unlist(.y), twosided = twosided))) %>%
                  tidyr::unnest(res) %>%
                  dplyr::select(-c(H0, H1)) %>%
                  dplyr::ungroup()

  list(roc = roc_results, raw = result, elapsed_time = elapsed_time)
}

#' Plot ROC curves
#' @description Plot ROC curves given results from \code{\link{simulate_roc}}
#' @param roc_results data frame \code{roc} from the output list of \code{\link{simulate_roc}}
#' @param method_names
#' @keywords plotting
#' @export
plot_roc <- function(roc_results, method_names=NULL) {
  if (!is.null(method_names)){
    roc_results = roc_results %>%
      dplyr::mutate(method = factor(method, levels = method_names))
  } else {
    roc_results = roc_results %>%
      dplyr::mutate(method = factor(method))
  }
  # customize line color
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p = roc_results %>%
    ggplot2::ggplot(ggplot2::aes(x = FPR, y = TPR, color = method)) +
    ggplot2::scale_colour_manual(values=cbbPalette)+
    ggplot2::geom_line(size=2)  +
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
