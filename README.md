### squidSim: a flexible simulation tool for linear mixed models.

<!-- [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/squid)](https://cran.r-project.org/package=squid)
[![Build
Status](https://travis-ci.org/squid-group/squid.svg?branch=master)](https://travis-ci.org/squid-group/squid)
[![Downloads](http://cranlogs.r-pkg.org/badges/squid?color=brightgreen)](https://cran.rstudio.com/package=squid)
[![total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/squid)](http://cranlogs.r-pkg.org/badges/grand-total/squid) -->

<img id='logo' src='./man/figures/squidSim_logo.png' align='left' alt='' style='padding-right:20px;' width='150'>

The `squidSim` package is a simulation tool that can be used for both research and educational purposes. Its main aim is to aid empiricists in the understanding of linear mixed effects model and hierarchical data. The package allows hierarchical data structure to be created or imported, which is then used as the basis for simulations. It also provides tools for sampling the data in multiple ways. The package has been designed with flexibility in mind and can simulate some a huge variety of different models, including the simulation of genetic and phylogenetic effects. 

The framework is suitable for performing simulation studies, determining optimal sampling designs for user-specific biological problems, and making simulation based inferences to aid in the interpretation of empirical studies. `squidSim` is also a teaching tool for biologists interested in learning, or teaching others, how to implement and interpret mixed-effects models.
<!-- Second, `squid` offers research opportunities to those who are already familiar with mixed-effects models, as `squid` enables the generation of datasets that users may download and use for a range of simulation-based statistical analyses such as power and sensitivity analysis of multilevel and multivariate data.
 -->


#### Installation

Currently the squidSim package is not on CRAN, but you can install the development version from GitHub using the devtools pacakge:

    # install.packages("devtools")
    devtools::install_github("squidgroup/squidSim")
    library(squidSim)


#### Getting started

We have an extensive [vignette](http://squidgroup.org/squidSim_vignette/), which we strongly recommend reading through before using the package.


#### Background to the project

<img id='logo' src='./man/figures/logo_2.png' align='left' alt='' width='160'>

**SQuID** stands for **S**tatistical **Qu**antification of **I**ndividual **D**ifferences and refers to a working group of behavioural and evolutionary ecologists (link to main group page).

Understanding the sources of phenotypic variance is a core component to a huge amount of research in ecology and evolutionary biology. Phenotypic variation is typically structured in a hierarchical way and the hierarchical modelling in mixed effect models provides a great tool to analyze and decompose such variation. Mixed models are very flexible statistical tools that provide a way to estimate the variation at these different levels, and represent the general statistical framework for ecology and evolutionary biology. Because of the progress in computational capacities mixed models have become increasingly popular among ecologists and evolutionary biologists over the last decade. However, fitting mixed model is not a straightforward exercise, and the way data are sampled among and within groups can have strong implications on the outcome of the model. Mixed models therefore present a highly flexible yet challenging tool.

A major aim of the SQuID working group is to provide educational resources to help empiricists learn more about the properties, potentials and interpretational challenges of mixed effects models. The group formed during a workshop 2013 with the aim of helping empiricists interested in decomposing phenotypic variance to get more familiar with the concept of hierarchical organization of traits, with mixed models and to avoid pitfalls caused by inappropriate sampling.

The first `squid` R package is fairly limited in the range of models that can be used to simulate data. The aim of the current package is to overcome these limitations to provide a coherent and reproducible framework with which to do simulations, both for research and education. 



