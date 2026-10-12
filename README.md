# FVDDPpkg


### R implementation of Fleming-Viot-dependent Dirchlet Processes

:computer: The package contains the first implementation of Fleming-Viot dependent Dirichlet processes, a model of Bayesian Nonparametric statistics used in populuation genetics.\
The model consists of a Hidden Markov Model (HMM) in which the latent signal is a distribution-valued stochastic process that takes the form of a finite mixture of Dirichlet processes, indexed by vectors that count how many times each value is observed in the population.

In particular, the package allows to:
* implement a Dirichlet process specifying
    * the baseline prior distribution $P_0$, which is described by its p.m.f (if it is discrete) or by its p.d.f. (it is is absolutely continuous), and by a function to sample from.
    * the intensity $\theta > 0$ associated to $P_0$.
* update the structure of the model with data, modifying the latent signal.
* propagate the latent signal in the past or in the future for any time, given the current state.
* smooth the latent signal given observations from the past, the present and the future.
* sample from the process.
* exploit the predictive structure to make inference.

For clarifications and examples of use, please refer to the vignette in the package documentation or on [GitHub](https://github.com/StefanoDamato/FVDDPpkg/blob/master/vignettes/FVDDPpkg.Rmd).

The package is described in [Ascolani, Damato and Ruggiero (2024)](https://doi.org/10.1089/cmb.2024.0600). For a rigorous theoretical approach, read [Ascolani et al. (2021)](https://projecteuclid.org/journals/bayesian-analysis/advance-publication/Predictive-inference-with-FlemingViot-driven-dependent-Dirichlet-processes/10.1214/20-BA1206.full) and [Ascolani et al. (2023)](https://projecteuclid.org/journals/bernoulli/volume-29/issue-2/Smoothing-distributions-for-conditional-FlemingViot-and-DawsonWatanabe-diffusions/10.3150/22-BEJ1504.short)


### Installation


:space_invader: The package has been released on CRAN!

To download the latest version, just run
```
install.packages("FVDDPpkg", dependencies = True)
```

:cat: To install from GitHub, make sure you have the `devtools` package installed.
```
install.packages("devtools")
```
and proceed to download using
```
devtools::install_github("StefanoDamato/FVDDPpkg", build_vignettes = TRUE)
```


### Citation

:book: If you use the package, please cite [Ascolani, Damato and Ruggiero (2024)](https://doi.org/10.1089/cmb.2024.0600), which presents the package and the algorithms it implements:
```
@article{AscolaniDamatoRuggiero2024,
author = {Filippo Ascolani and Stefano Damato and Matteo Ruggiero},
title = {{An R Package for Nonparametric Inference on Dynamic Populations with Infinitely Many Types}},
volume = {31},
journal = {Journal of Computational Biology},
number = {12},
pages = {1305 -- 1311},
year = {2024},
doi = {10.1089/cmb.2024.0600},
URL = {https://doi.org/10.1089/cmb.2024.0600}
}
```


### Contacts

:e-mail: For questions, contributions or bug reports contact me at `stefano.damato@idsia.ch` or use the [issues section](https://github.com/StefanoDamato/FVDDPpkg/issues) on the GitHub repository.
