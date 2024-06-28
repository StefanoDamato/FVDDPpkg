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

For a rigorous theoretical approach, read [Ascolani et al. (2021)](https://projecteuclid.org/journals/bayesian-analysis/advance-publication/Predictive-inference-with-FlemingViot-driven-dependent-Dirichlet-processes/10.1214/20-BA1206.full) and [Ascolani et al. (2023)](https://projecteuclid.org/journals/bernoulli/volume-29/issue-2/Smoothing-distributions-for-conditional-FlemingViot-and-DawsonWatanabe-diffusions/10.3150/22-BEJ1504.short)


### Installation


:cat: To install from GitHub, make sure you have the `devtools` package installed.
```
install.packages("devtools")
```
and proceed to download using
```
devtools::install_github("StefanoDamato/FVDDPpkg", build_vignettes = TRUE)
```

:space_invader: Soon too appear on CRAN!

After the acceptance, just run
```
install.packages("FVDDPpkg", dependencies = True)
```

### Citation

:book: A paper on the package and new algorithms it implements is currently under review. Its publication will be notified here now.

At the moment, the most comprehensive works you can cite are
```
@article{AscolaniLijoiRuggiero2021,
author = {Filippo Ascolani and Antonio Lijoi and Matteo Ruggiero},
title = {{Predictive inference with Fleming–Viot-driven dependent Dirichlet processes}},
volume = {16},
journal = {Bayesian Analysis},
number = {2},
publisher = {International Society for Bayesian Analysis},
pages = {371 -- 395},
keywords = {Chinese restaurant, conveyor belt, generalized Pólya urn, Hidden Markov model, predictive distribution, random partition},
year = {2021},
doi = {10.1214/20-BA1206},
URL = {https://doi.org/10.1214/20-BA1206}}
```
and
```
@article{AscolaniLijoiRuggiero2023,
author = {Filippo Ascolani and Antonio Lijoi and Matteo Ruggiero},
title = {{Smoothing distributions for conditional Fleming–Viot and Dawson–Watanabe diffusions}},
volume = {29},
journal = {Bernoulli},
number = {2},
publisher = {Bernoulli Society for Mathematical Statistics and Probability},
pages = {1410 -- 1434},
keywords = {Dirichlet process, Duality, gamma random measures, Hidden Markov model, optimal filtering, prediction},
year = {2023},
doi = {10.3150/22-BEJ1504},
URL = {https://doi.org/10.3150/22-BEJ1504}
}
```


### Contacts

:e-mail: For questions, contributions or bug reports contact me at `stefano.damato@idsia.ch` or use the [issues section](https://github.com/StefanoDamato/FVDDPpkg/issues) on the GitHub repository.