# nett 


`nett` is an R package for the analysis of network data with a focus on **community detection** and associated models.

Check out the articles for some examples of how to use the package.

The package is under development and hopefully more features will be added as time goes by.

Some of the implemented functionality are follows:
- Spectral clustering with regularization
- Conditional pseudo-likelihood for community detection ([Amini, Chen, Bickel and Levina](https://projecteuclid.org/euclid.aos/1382547514)).
- Spectral goodness-of-test for SBM and DCSBM (inspired by Bickel and Sarkar, and [Lei](https://projecteuclid.org/euclid.aos/1452004791)'s work). 
- Likelihood ratio tests and BIC selection for SBM and DCSBM (inspired by [Bickel and Wang](https://projecteuclid.org/euclid.aos/1494921948)'s work among others.)
- Likelihood computations for SBM and DCSBM.
- Network Adjusted Chi-square test for SBM and DCSBM ([Zhang and Amini](https://arxiv.org/abs/2012.15047)).
- Bethe-Hessian Selection for DCSBM (inspired by [Le and Levina](https://arxiv.org/abs/1507.00827)'s work).
- ...

Most of the computations haven been adapted to Poisson models to make them fast and scalable. 

## Installation
To install, you can use the following command

<pre>
devtools::install_github("aaamini/nett")
</pre>

To also build the vignettes when installing use the follwing:

<pre>
devtools::install_github("aaamini/nett", build_vignettes = TRUE)
</pre>

This will take some time. Alternatively, just refer to the articles on this website (which are the the built vignettes).

### Related 
See the related repository [linfanz/nac-test](https://github.com/linfanz/nac-test), for some experiments comparing goodness-of-fit and model selection approaches.
