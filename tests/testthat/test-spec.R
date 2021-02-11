n = 1500
Ktru = 4
lambda = 15 # expected average degree
oir = 0.1
pri = 1:Ktru

set.seed(1234)
theta <- EnvStats::rpareto(n, 2/3, 3)
B = pp_conn(n, oir, lambda, pri=pri, theta)$B
z = sample(Ktru, n, replace=T, prob=pri) # randomly smaple "true community labels"
A = sample_dcsbm(z, B, theta) # sample the adjacency matrix
zh = spec_clust(A, K=4)

test_that("spec works -- pp", {
  expect_equal(signif(compute_mutual_info(z, zh),5), 0.84595)
})

set.seed(1400)
theta <- EnvStats::rpareto(n, 2/3, 3)
B = gen_rand_conn(n, Ktru, lambda = lambda, gamma = 0.1, pri=pri)
z = sample(Ktru, n, replace=T, prob=pri) # randomly smaple "true community labels"
A = sample_dcsbm(z, B, theta) # sample the adjacency matrix

test_that("spec works -- gen_rand", {
  expect_equal(signif(compute_mutual_info(z, spec_clust(A, K=4)),5), 0.92094)
  expect_equal(signif(compute_mutual_info(z, spec_clust(A, K=4, tau=0.1)),5), 0.92781)
})
