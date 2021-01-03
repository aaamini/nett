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
#' \item{roc}{}
#' \item{raw}{}
#' \item{elapsed_time}{symstem elapsed time for generating ROC data}
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
    #   tibble::enframe(apply_methods(null_data), name="method", value = "tstat") %>% # turns a named list into a two-column tibble
    #     tidyr::unnest(tstat) %>%
    #     tibble::add_column(H=0),
    #   tibble::enframe(apply_methods(alt_data), name="method", value = "tstat") %>%
    #     tidyr::unnest(tstat) %>%
    #     tibble::add_column(H=1)
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

  # method_names = unique(result$method)
  # n_methods = length(method_names)
  # roc_results <- NULL
  # for (r in 1:n_methods) {
  #   mtd_name <- method_names[r]
  #   curr_result <- dplyr::filter(result, method == mtd_name)
  #   T0 <- dplyr::filter(curr_result, H==0)$tstat
  #   twosided <- curr_result$twosided[1] # all elements should be the same
  #   T1 <- dplyr::filter(curr_result, H==1)$tstat
  #
  #   roc_results <- rbind(roc_results,
  #                        get_roc(T0, T1, twosided = twosided) %>%
  #                          tibble::add_column(method=mtd_name) )
  # }
  list(roc = roc_results, raw = result, elapsed_time = elapsed_time)
}

# simulate_roc = function(methods, gen_null_data, gen_alt_data, nruns = 100,
#                         twosided = T,
#                         core_count = parallel::detectCores() - 1) {
#   # method_names <- gsub("p","+",names(methods))
#   method_names <- names(methods)
#   n_methods <- length(methods)
#   if (length(twosided) == 1) {
#     twosided = rep(twosided, n_methods)
#   }
#
#   simulate_roc_run = function(j) {
#     null_data = gen_null_data()
#     alt_data = gen_alt_data()
#
#     result = NULL
#     for (r in 1:n_methods) {
#       mtd_name = method_names[r]
#       mtd_fun = methods[[r]]
#       result = dplyr::bind_rows(result, tibble::tibble(method = mtd_name, H=0, tstat=mtd_fun(null_data)) )
#       result = dplyr::bind_rows(result, tibble::tibble(method = mtd_name, H=1, tstat=mtd_fun(alt_data)) )
#     }
#     result
#   }
#
#   elapsed_time <- as.numeric(
#     system.time(
#       result <- do.call(rbind, parallel::mclapply(1:nruns, simulate_roc_run, mc.cores = core_count))
#     )["elapsed"]
#   )
#
#   roc_results <- NULL
#   for (r in 1:n_methods) {
#     mtd_name <- method_names[r]
#     curr_result <- dplyr::filter(result, method == mtd_name)
#     T0 <- dplyr::filter(curr_result, H==0)$tstat
#     T1 <- dplyr::filter(curr_result, H==1)$tstat
#
#     roc_results <- rbind(roc_results,
#                          get_roc(T0, T1, twosided = twosided) %>%
#                            tibble::add_column(method=mtd_name) )
#   }
#   list(roc=roc_results, raw=result, elapsed_time=elapsed_time)
# }

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
