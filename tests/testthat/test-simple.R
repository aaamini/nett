library(Matrix)
library(EnvStats)

n = 1000
Ktru = 4
lambda = 10
oir = 0.1

set.seed(1)
theta <- rpareto(n, 2/3, 3)
temp = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)
A = temp$adj
z = temp$labels
B = temp$B

test_that("spectral clustering works", {
  zh = spec_clust(A,4)
  expect_equal(round(compute_mutual_info(z, zh),2), 0.59)
})

test_that("spectral test", {
  zh = spec_clust(A,4)
  expect_equal(round(spec_test(A, 4),2), 121.9)
})

test_that("estim_dcsbm", {
  # out = estimDCSBM(A, z)
  out = estim_dcsbm(A, z)
  expect_equal(round(norm(B - out$B) / norm(B), 3), 0.247)
  expect_equal(round(max(abs(out$theta-theta) / theta),3), 2.214)
})

# library(microbenchmark)
# n = 2000
# Ktru = 3
# lambda = 10
# theta <- rpareto(n, 2/3, 3)
# temp = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)
# A = temp$adj
# z = temp$labels
# B = temp$B
#
# microbenchmark(forloop = computeBlockSums2(A,z), Matrix = computeBlockSums(A,z), setup = {A = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)$adj});
# microbenchmark(outer = outer(ns,ns), mprod = ns %*% t(ns), setup = {ns <- as.vector(table(z))})
#
# microbenchmark("estim_sbm" = estim_sbm(A,z), "estimSBM" =  estimSBM(A, z),
#                setup = {A = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)$adj}, times = 20)
#
# microbenchmark("estim_dcsbm" = estim_dcsbm(A,z), "estimDCSBM" =  estimDCSBM(A, z),
#                setup = {theta <- rpareto(n, 2/3, 3); A = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)$adj}, times = 20)
#
#
# ns <- as.vector(table(z))
#
# norm(B - estimSBM(A, z)) / norm(B)
# norm(B - estim_sbm(A, z)) / norm(B)
#
# computeBlockSums2(A,z)
# sum(A[z==1,z==3])
#
# #Ktru=3
# zh = spec_clust(A,Ktru)
# estimSBM(A, zh)
# estim_sbm(A, zh)
#
#
# out = estimDCSBM(A, zh)
# out$B
# out$theta
#
# out2 = estim_dcsbm(A, zh)
# out2$B
# out2$theta
#
# out$B - out2$B
# out$theta - out2$theta
#
# # l2_norm = function(x) sqrt(sum(x^2))
# round(max(abs((out$theta-theta) / theta)),2)
# hist(abs(out$theta-theta) / theta)
#
