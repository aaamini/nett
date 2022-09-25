#' @import stats
customSampleFun <- function(m, n) {
  # sample n numbers from 1 : m
  # sample() is too slow when m is large
  while (1) {
    x = ceiling(runif(n*2) * m)
    y = unique(x)
    if (length(y) >= n) {
      break
    }
  }
  return(sample(y, n))
}

# MySampleFun <- customSampleFun
# MySampleFun <- function(m,n) sample.int(m, n, useHash = TRUE) # gives error when n > m/2
MySampleFun <- function(m,n) igraph::sample_seq(1,m,n)

# check MySampleFun
# MySampleFun <- sample
