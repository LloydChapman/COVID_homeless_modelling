# Comparison of interventions against COVID-19 outbreaks in homeless shelters

This repository contains Approximate Bayesian Computation (ABC) code and simulation code for the analyses in 'Comparison of infection control strategies to reduce COVID-19 outbreaks in homeless shelters in the United States: a simulation study'. The code implements a discrete-time stochastic SEIR simulation model of COVID-19 transmission in a closed environment (here a homeless shelter) with importation of infection from the local community. The model is fitted to data on numbers of PCR-positive and negative individuals  from outbreaks in 5 homeless shelters in San Francisco, Boston and Seattle, and used to predict the impact of different intervention strategies on the probability of averting an outbreak over 30 days in a representative homeless shelter into which a single latently infected individual is introduced.

## Prerequisites

* [R version 4.0.0](https://www.r-project.org/)

* The following R packages are required to run the code:
  * mvnfast
  * actuar
  * ggplot2
  * reshape2
  * Hmisc
  * gsubfn
  * gridExtra
  * doParallel
  * abind

## Data

All data required to run the code is available in the [data](data) subfolder.

## Installing

Clone/download this project into a folder on your machine using the green button at the top right of this page.

### Running the code

The required R packages can be installed by running the following line of code in R

```R
> install.packages(c("mvnfast","actuar","ggplot2","reshape2","Hmisc","gsubfn","gridExtra","doParallel","abind"))
```

The model calibration can then be run in R by entering

```R
> source("run_calibration.R")
```

at the command prompt, or by navigating to the downloaded code folder in a terminal window on Mac/Linux and entering

```
% Rscript run_calibration.R
```
 
or in Windows command line by entering

```
C:\>"C:\<path>\<to>\Rscript.exe" C:<path>\<to>\run_calibration.R
```

The intervention simulations and sensitivity analysis can be run in R similarly with

```R
> source("run_interventions.R")
> source("run_sensitivity_analysis.R")
```

or via the command line (Mac/Linux) with

```
% Rscript run_interventions.R
% Rscript run_interventions.R
```

or from the Windows command line with

```
C:\>"C:\<path>\<to>\Rscript.exe" <path>\<to>\run_interventions.R
C:\>"C:\<path>\<to>\Rscript.exe" <path>\<to>\run_interventions.R
```

The SEIR transmission model is implemented in the COVID_homeless_functions.R script, and the fixed model parameters are set in set_nat_hist_pars.R. The ABC Sequential Monte Carlo (SMC) algorithm can be found in ABC_SMC.R.

## Built With

* [R version 4.0.0 (2020-04-24)](https://www.r-project.org/)

## Author

* Lloyd Chapman: <lloyd.chapman@ucsf.edu>

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments

The ABC SMC code adapts the code for case study 2 in [1] available [here](https://github.com/amanda-minter/abc_R/tree/master/case_2) to enable fitting to both discrete and continuous parameters following the ABC model-selection algorithm described in [2,3].

## References
1. Minter A, Retkute R. Approximate Bayesian Computation for infectious disease modelling. Epidemics. 2019;29:100368. doi:[10.1016/j.epidem.2019.100368](https://doi.org/10.1016/j.epidem.2019.100368)

2. Toni T, Welch D, Strelkowa N, Ipsen A, Stumpf MPH. Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems. J R Soc Interface. 2009;6(31):187–202. doi:[10.1098/rsif.2008.0172](https://royalsocietypublishing.org/doi/10.1098/rsif.2008.0172)

3. Toni T, Stumpf MPH. Simulation-based model selection for dynamical systems in systems and population biology. Bioinformatics. 2009;26(1):104–10. doi:[10.1093/bioinformatics/btp619](https://doi.org/10.1093/bioinformatics/btp619)
