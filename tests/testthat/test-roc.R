

methods = list(ttest = function(x) mean(x),
               wilcox = function(x) wilcox.test(x)$statistic  )

n = 10
# set.seed(1)
out =  simulate_roc(methods,
                    gen_null_data = function() rnorm(n),
                    gen_alt_data = function() rnorm(n, mean=0.5),
                    nruns = 20)

test_that("roc works", {
  expect_equal(dim(out$roc)[2], 3)
})

# plot_roc(out$roc)
# table(out$roc %>% dplyr::filter(method == "wilcox") %>% dplyr::select(FPR))
