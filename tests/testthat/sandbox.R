library(Matrix)
library(EnvStats)

n = 1000
Ktru = 4
lambda = 10
oir = 0.1

set.seed(1)
theta <- EnvStats::rpareto(n, 2/3, 3)
#n*mean(theta/max(theta))^2
temp = quickDCSBM(n, lambda,  Ktru, oir = oir, theta, pri=1:Ktru, normalize_theta = F)
A = temp$adj
z = temp$labels
B = temp$B
#theta = temp$theta
mean(rowSums(A))

# library(igraph)
# polblogs2 = read_graph("../CAC test/code/new R code/data/polblogs2.gml", format = "gml")
# usethis::use_data(polblogs)

# gen_rand_conn = function(n, K, lambda, gamma, pri = rep(1,K)/K) {
#   B = matrix(runif(K^2),K)
#   B = (B + t(B))/2
#   rand_perm_mat = diag(K)[, sample(K)]
#   B = (1-gamma)*rand_perm_mat + gamma*B
#   scale = get_dcsbm_exav_deg(n, pri, B)
#   B*lambda/scale
# }
#
# B = gen_rand_conn(n, Ktru, lambda, gamma)
# A <- fastDCSBM(z, B, theta = theta)
# compute_mutual_info(z, spec_clust(A, Ktru))

# fastDCSBM(z, B, theta)
# Z = label_vec2mat(z)
# P = diag(theta) %*% Z %*% B %*% t(Z) %*% diag(theta)
# A = Matrix(matrix(runif(n^2),n) < P)
# mean(rowSums(A))

