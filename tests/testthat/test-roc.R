

# methods = list(ttest = function(x) mean(x),
#                wilcox = function(x) wilcox.test(x)$statistic  )

apply_methods = function(x) {
  tibble::tribble(
    ~method, ~tstat, ~twosided,
    "ttest", as.numeric(t.test(x)$statistic), T, #mean(x), T,
    "wilcox", as.numeric(wilcox.test(x)$statistic), T
  )
}

n = 10
out =  simulate_roc(apply_methods,
                    gen_null_data = function() rnorm(n),
                    gen_alt_data = function() rnorm(n, mean=0.5),
                    core_count = 1,
                    nruns = 500, seed = 15)
printf('time = %3.3f', out$elapsed_time)
plot_roc(out$roc)

test_that("roc works", {
  expect_equal(dim(out$roc)[2], 4)
})


# plot_roc(out$roc)
# table(out$roc %>% dplyr::filter(method == "wilcox") %>% dplyr::select(FPR))


# K <- 4 # number of communities
# n <- 500
# oir <- 0.05 # out-in ratio
# lambda <- 10 # average degrees
# theta <- EnvStats::rpareto(n, 2/3, 3)
#
# gen_null_data = function() quickDCSBM(n, lambda,  K_H0, oir = oir, theta)$adj
# gen_alt_data = function() quickDCSBM(n, lambda,  K_H1, oir = oir, theta)$adj
#
# K_H0 <- K
# K_H1 <- K-1
# apply_methods = function(A) {
#   z0 = spec_clust(A, K_H0)
#   z0p = spec_clust(A, K_H0+1)
#   z1 = spec_clust(A, K_H1)
#   tibble::tribble(
#     ~method, ~tstat, ~twosided,
#     "NAC+", CCTest(A, K_H0, z=z0, y=z0p)$total, T,
#     "SNAC+", CCsub(A, K_H0, z=z0, ratio=.5, fromEachCommunity=TRUE)$total, T,
#     "DCSC", spec_test(A, K_H0, z=z0), T,
#     "LR", eval_dcsbm_loglr(A, cbind(z0, z1), poi=T), F
#   )
# }
#
# apply_methods(gen_null_data())
# out = simulate_roc(apply_methods,
#                    gen_null_data,
#                    gen_alt_data,
#                    nruns = 20, core_count = 1, seed = 1)
# printf('time = %3.3f', out$elapsed_time)
# plot_roc(out$roc)

# temp = out$raw %>% dplyr::filter(method == "SNAC+")
# temp %>% ggplot2::ggplot(ggplot2::aes(x=tstat, color=as.factor(H))) + ggplot2::geom_histogram()
