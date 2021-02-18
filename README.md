# BayesOnline_ChangePoint_Detection
This is an implementation of [Bayesian Online ChangePoint Detection](https://arxiv.org/abs/0710.3742) as described in the paper and is further extended for finding change-points in an AR process.

The Q1.m file contains the implementation of bocd for well-drilling nuclear magnetic
response data (similar to the one used in the paper). 
Q2.m employs bocd and RLS in tandem to predict the changepoints of the variance of the
signal's innovations.

* Prior: Normal-Gamma
* Likelihood: Gaussian
* Hazard Function: Geometric Distribution

Most implementation details are in report.pdf

## Some Useful links for understanding the paper:
1.	http://gregorygundersen.com/blog/2019/08/13/bocd/
2.	http://gregorygundersen.com/blog/2020/10/20/implementing-bocd/
3.  [Murphy- to evaluate the posterior predictive](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf)
