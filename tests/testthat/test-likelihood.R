# library(Matrix)
# library(EnvStats)
#
# n = 400
# Ktru = 4
# lambda = 10
# oir = 0.1
#
# set.seed(1)
# theta <- rpareto(n, 2/3, 3)
# temp = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)
# A = temp$adj
# z = temp$labels
# B = temp$B
#
# test_that("likelihood computation", {
#   expect_equal( round(eval_dcsbm_like(A, z, poi = F), 3), -6761.271 )
#   expect_equal( round(eval_dcsbm_like(A, z, poi = T), 3), -6802.076 )
# })
#
# # test_that("likelihood ratio", {
# #   z0 = spec_clust(A,4)
# #   z1 = spec_clust(A,5)
# #   expect_equal( round(eval_dcsbm_loglr(A, cbind(z0,z1), poi=F),5), 1.00471)
# #   expect_equal( round(eval_dcsbm_loglr(A, cbind(z0,z1), poi=T),5), 1.00513)
# # })
#
#
#
# # n <- nrow(A)
# #
# # d <- colSums(A)
# # upper.index <- which(upper.tri(A))
# # a <- A[upper.index]
# #
# # K = Ktru
# # g <- z
# # n.K <- as.numeric(table(g))
# # Pi <- n.K/n
# # B.star <- matrix(0,K,K)
# #
# # for(j in 1:(K)){
# #   for(k in j:K){
# #     j.index <- which(g==j)
# #     k.index <- which(g==k)
# #     B.star[j,k] <- sum(A[j.index,k.index])
# #     B.star[k,j] <- B.star[j,k]
# #   }
# # }
# # b.row <- rowSums(B.star)
# # Z <- matrix(0,n,K)
# # Z[cbind(1:n,g)] <- 1
# # Z.b.row <- Z%*%b.row
# # theta <- d/Z.b.row
# # Z.theta <- Z*as.numeric(theta)
# # dc.P <- Z.theta%*%B.star%*%t(Z.theta)
# # dc.p <- dc.P[upper.index]
# # dc.p[dc.p>(1-1e-6)] <- 1-1e-6
# # dc.p[dc.p<1e-6] <- 1e-6
# # (dcsbm.ll <- sum(a*log(dc.p)) + sum((1-a)*log(1-dc.p)) + sum(n.K*log(Pi)))
# #
# # eval_dcsbm_like(A, z, poi = F)
#
