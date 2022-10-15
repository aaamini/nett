# nett 

`nett` is an R package for the analysis of network data with a focus on **community detection** and implements multiple methods for hypothesis testing. It includes model selection and goodness-of-fit tests of SBM/DCSBM to network data, which are useful in network statistical analysis. Some of the implemented functionality are follows:
- Spectral clustering with regularization
- Conditional pseudo-likelihood for community detection ([Amini, Chen, Bickel and Levina](https://projecteuclid.org/euclid.aos/1382547514)).
- Spectral goodness-of-test for SBM and DCSBM (inspired by [Bickel and Sarkar](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12117), and [Lei](https://projecteuclid.org/euclid.aos/1452004791)'s work). 
- Likelihood ratio tests and BIC selection for SBM and DCSBM (inspired by [Wang and Bickel](https://projecteuclid.org/euclid.aos/1494921948)'s work among others.)
- Likelihood computations for SBM and DCSBM.
- Network Adjusted Chi-square test for SBM and DCSBM ([Zhang and Amini](https://arxiv.org/abs/2012.15047)).
- Bethe-Hessian Selection for DCSBM (inspired by [Le and Levina](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-16/issue-1/Estimating-the-number-of-communities-by-spectral-methods/10.1214/21-EJS1971.full)'s work).
- ...

Most of the computations haven been adapted to Poisson models to make them fast and scalable. 

Check out the [articles](https://aaamini.github.io/nett/) for some examples of how to use the package.

## Installation
To install, you can use the following command

<pre>
devtools::install_github("aaamini/nett")
</pre>

## Related repo
See the related repo [linfanz/nac-test](https://github.com/linfanz/nac-test), for some experiments comparing goodness-of-fit and model selection approaches.
