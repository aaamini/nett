### declare global variable
globalVariables(c("K", "value", "fxx"))


#' Resampled SNAC+
#'
#' Compute SNAC+ with resampling
#' @param A adjacency matrix
#' @param nrep number of times SNAC+ is computed
#' @param Kmin minimum community number to use in SNAC+
#' @param Kmax maximum community number to use in SNAC+
#' @param ncores number of cores to use in the parallel computing
#' @param seed seed for random sampling
#' @return A data frame with columns specifying repetition cycles,
#' number of community numbers and the value of SNAC+ statistics
#' @export
snac_resample = function(A, nrep = 20, Kmin = 1, Kmax = 13,
                         ncores =  parallel::detectCores()-1,
                         seed = 1234) {
    Ks = Kmin:Kmax

    # parallel version does not work under windows
    # cl <- parallel::makeForkCluster(ncores)
    # cl <- parallel::makeCluster(ncores)
    # doParallel::registerDoParallel(cl)
    # doRNG::registerDoRNG(seed)
    # labels = sapply(Kmin:(Kmax+1), function(k) spec_clust(A, k))
    # Tstat = do.call(rbind,
    #     foreach::foreach(t = 1:nrep) %dopar% {
    #         data.frame(
    #             itr = rep(t, length(Ks)),
    #             K = Ks,
    #             value = sapply(Ks, function(k) snac_test(A, k, labels[ , k-Kmin+1])$stat)
    #         )
    #     })
    #
    # # Tstat = bind_rows(Tstat)
    # parallel::stopCluster(cl)

    labels = sapply(Kmin:(Kmax+1), function(k) spec_clust(A, k))
    Tstat = do.call(rbind, parallel::mclapply(1:nrep,
                                              function(t){
                                                data.frame(
                                                  itr = rep(t, length(Ks)),
                                                  K = Ks,
                                                  value = sapply(Ks, function(k) snac_test(A, k, labels[ , k-Kmin+1])$stat)
                                                )
                                              },
                                              mc.cores = ncores
                                              )
                      )

    # Tstat = bind_rows(Tstat)
    # parallel::stopCluster(cl)


    Tstat
}

# fp = f', fpp = f''
detect_dip = function(fp){
    which(fp > 0)[1]
}

# compute_curvature = function(fp, fpp) {
#     abs(fpp) / ((1 + fp^2)^(1.5))
# }

detect_elbow = function(curvature) {
    which.max(curvature)[1]
}

fit_ss = function(x, y, xx, spar = NULL, trunc_type = "none") {
    trunc_fun = switch(trunc_type, "none" = function(x) signif(x, 2),
                       "floor" = floor, "ceil" = ceiling, "round" = round)

    fitted_ss = smooth.spline(x, y, spar = spar)
    fxx_ss = predict(fitted_ss, xx, deriv = 0)$y
    fp_ss = predict(fitted_ss, xx, deriv = 1)$y
    fpp_ss = predict(fitted_ss, xx, deriv = 2)$y

    list(dip = trunc_fun( xx[detect_dip(fp_ss)]),
         elbow1 = trunc_fun( xx[detect_elbow(fpp_ss)] ),
         fxx = fxx_ss)
}

#' Plot community profiles
#'
#' Plot the smooth community profiles based on a resampled statistic
#' @param tstat dataframe that has a column 'value' as statistic in the plot
#' and a column 'K' as its corresponding community number
#' @param net_name name of network
#' @param trunc_type method to round the dip/elbow point as the estimated
#' community number
#' @param  spar the sparsity level of fitting spline to the value of \code{tstat}
#' @param plot_null_spar whether to plot the spline with zero sparsity
#' @param alpha transparency of the points in the plot
#' @param base_font_size font size of the plot
#' @return smooth profile plot of a network
#' @export
plot_smooth_profile = function(
    tstat, net_name = "", trunc_type = "none",
    spar=0.3, plot_null_spar = TRUE,
    alpha = 0.3, # transparency
    base_font_size = 12
) {
    trunc_fun = switch(trunc_type, "none" = function(x) signif(x, 2),
                       "floor" = floor, "ceil" = ceiling, "round" = round)

    y = tstat$value
    x = tstat$K
    Kmax = max(x)
    xx <- seq(1, Kmax, length.out = 1000)

    ss_res1 = fit_ss(x, y, xx, spar = spar, trunc_type = trunc_type)
    color2 = "#FB674B"

    optK_label = c(shQuote(net_name), bquote(K %~~% .(ss_res1$elbow1) ~","~.(ss_res1$dip)))
    p = ggplot2::ggplot(tstat, ggplot2::aes(x=K, y=value)) +
        ggplot2::geom_point(size = 2, color="#6F98C2", alpha = alpha) +
        ggplot2::geom_line(ggplot2::aes(x=xx, y=fxx), data = data.frame(xx=xx, fxx=ss_res1$fxx), size = 1.5, color="black")  +
        ggplot2::scale_x_continuous(breaks = unique(tstat$K)) +
        ggplot2::ylab("SNAC+") +
        ggplot2::theme_bw(base_size = base_font_size)
    if (!is.null(spar) && plot_null_spar) {
        ss_res2 = fit_ss(x, y, xx, spar = NULL, trunc_type = trunc_type)
        p = p +
            ggplot2::geom_line(ggplot2::aes(x=xx, y=fxx), data = data.frame(xx=xx, fxx=ss_res2$fxx), size = 1.5, color = color2, linetype = "dashed")
        optK_label =c(optK_label, bquote(K %~~% .(ss_res2$elbow1) ~","~.(ss_res2$dip)))
    }
    p + ggplot2::annotate("text", x=Kmax-1.75, y=max(tstat$value)*c(0.95,0.9,0.85),
                 label = optK_label, size = 5, parse =TRUE, color = c("black", "black", color2))
}





