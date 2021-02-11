# library(Matrix)
# library(EnvStats)
#
# n = 2000
# Ktru = 4
# lambda = 20
# oir = 0.1
#
# set.seed(1)
# pri = 1:Ktru
# theta <- EnvStats::rpareto(n, 2/3, 3)
# summary(theta)
# B = pp_conn(n, oir, lambda, pri=pri, theta)$B
# z = sample(Ktru, n, replace=T, prob=pri)
# #n*mean(theta/max(theta))^2
# # temp = quickDCSBM(n, lambda,  Ktru, oir = oir, theta, pri=1:Ktru, normalize_theta = F)
# # A = temp$adj
# # z = temp$labels
# # B = temp$B
# #theta = temp$theta
# # mean(rowSums(A))
# A = sample_dcsbm(z, B, theta)
# d = rowSums(A)
# hist(d, 40)
# summary(d)
# mean(rowSums(A))
# Matrix::image(A)
# gr = igraph::graph_from_adjacency_matrix(A, "undirected")
# out = plot_net(gr, community = z)
#
# plot_deg_dist(gr, log=T)
# summary(igraph::degree(out$gr))
#
# plot_net(polblogs, community = igraph::V(polblogs)$community)
# plot_deg_dist(polblogs, logx=T)
#
# check_pkg_and_stop("ig")
#
# idx = order(theta, decreasing = T)
# head(cbind(rowSums(A2)[idx], theta[idx]), 20)
#
# tstat = snac_resample(A2, nrep = 10, ncores = 3)
# plot_smooth_profile(tstat, "temp", trunc_type = "none", spar=0.3, plot_null_spar = T)
#
#
# A = igraph::as_adj(igraph::as.undirected(polblogs))
# tstat = snac_resample(A, nrep = 30, ncores = 3)
# plot_smooth_profile(tstat, "temp", trunc_type = "none", spar=0.3, plot_null_spar = T)
#
#
# n = 1000
# Ktru = 4
# lambda = 20
# d = Ktru
# labels = sample(Ktru, n, replace=T, prob=rep(1,Ktru))
# labels = sort(labels)
# mu = diag(Ktru)
# z = 2*mu[labels, ] + 0.75*matrix(rnorm(n*d), n)
#
#
# plot(z[,1],z[,2])
#
# theta =  EnvStats::rpareto(n, 3/4, 4)
#
# # theta = rep(1,n)
# # cbind(pair_dist2_mat(z)[1,], labels)
#
# A = sample_dclvm(z, lambda, theta)
#
# Matrix::image(A)
# gr = igraph::graph_from_adjacency_matrix(A, "undirected")
# summary(igraph::degree(gr))
# plot_deg_dist(gr)
# plot_net(gr, extract_lcc = T, community = labels)
# igraph::components(gr)$csize
#
#
# # library(Matrix)
# # A = Matrix(matrix(runif(n^2),n) <  0.1*outer(theta,theta)*exp(-as.matrix(dist(z))^2))*1
# # A[upper.tri(A)] = 0
# # diag(A) = 0
# # A = (A + t(A))
# # image(A)
# # isSymmetric(A)
#
# # microbenchmark::microbenchmark(quickDCSBM(n, lambda,  Ktru, oir = oir, theta, pri=1:Ktru, normalize_theta = F),
# #                                sample_dcsbm(z, B, theta), times = 10)
#
# # gen_rand_conn = function(n, K, lambda, gamma, pri = rep(1,K)/K) {
# #   B = matrix(runif(K^2),K)
# #   B = (B + t(B))/2
# #   rand_perm_mat = diag(K)[, sample(K)]
# #   B = (1-gamma)*rand_perm_mat + gamma*B
# #   scale = get_dcsbm_exav_deg(n, pri, B)
# #   B*lambda/scale
# # }
# #
# # B = gen_rand_conn(n, Ktru, lambda, gamma)
# # A <- fastDCSBM(z, B, theta = theta)
# # compute_mutual_info(z, spec_clust(A, Ktru))
#
# # fastDCSBM(z, B, theta)
# # Z = label_vec2mat(z)
# # P = diag(theta) %*% Z %*% B %*% t(Z) %*% diag(theta)
# # A = Matrix(matrix(runif(n^2),n) < P)
# # mean(rowSums(A))
#
#
#
# test_that("spectral clustering works", {
#   zh = spec_clust(A, 4, ignore_first_col = F)
#   zh2 = fastCPL(A, 4, ilabels = zh)
#
#   expect_equal(round(compute_mutual_info(z, zh),2), 0.55)
#   expect_equal(round(compute_mutual_info(z, zh2),2), 0.62)
#
#   eval_dcsbm_bic(A, zh, 4)
# })
#
# test_that("spectral test", {
#   zh = spec_clust(A,4)
#   expect_equal(round(spec_test(A, 4),2), 121.9)
# })
#
# test_that("estim_dcsbm", {
#   # out = estimDCSBM(A, z)
#   out = estim_dcsbm(A, z)
#   expect_equal(round(norm(B - out$B) / norm(B), 3), 0.247)
#   expect_equal(round(max(abs(out$theta-theta) / theta),3), 2.214)
# })
#
# # library(microbenchmark)
# # n = 2000
# # Ktru = 3
# # lambda = 10
# # theta <- rpareto(n, 2/3, 3)
# # temp = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)
# # A = temp$adj
# # z = temp$labels
# # B = temp$B
# #
# # microbenchmark(forloop = computeBlockSums2(A,z), Matrix = computeBlockSums(A,z), setup = {A = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)$adj});
# # microbenchmark(outer = outer(ns,ns), mprod = ns %*% t(ns), setup = {ns <- as.vector(table(z))})
# #
# # microbenchmark("estim_sbm" = estim_sbm(A,z), "estimSBM" =  estimSBM(A, z),
# #                setup = {A = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)$adj}, times = 20)
# #
# # microbenchmark("estim_dcsbm" = estim_dcsbm(A,z), "estimDCSBM" =  estimDCSBM(A, z),
# #                setup = {theta <- rpareto(n, 2/3, 3); A = quickDCSBM(n, lambda,  Ktru, oir = oir, theta)$adj}, times = 20)
# #
# #
# # ns <- as.vector(table(z))
# #
# # norm(B - estimSBM(A, z)) / norm(B)
# # norm(B - estim_sbm(A, z)) / norm(B)
# #
# # computeBlockSums2(A,z)
# # sum(A[z==1,z==3])
# #
# # #Ktru=3
# # zh = spec_clust(A,Ktru)
# # estimSBM(A, zh)
# # estim_sbm(A, zh)
# #
# #
# # out = estimDCSBM(A, zh)
# # out$B
# # out$theta
# #
# # out2 = estim_dcsbm(A, zh)
# # out2$B
# # out2$theta
# #
# # out$B - out2$B
# # out$theta - out2$theta
# #
# # # l2_norm = function(x) sqrt(sum(x^2))
# # round(max(abs((out$theta-theta) / theta)),2)
# # hist(abs(out$theta-theta) / theta)
# #
