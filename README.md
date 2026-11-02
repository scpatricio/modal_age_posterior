[![arXiv shield](https://img.shields.io/badge/arXiv-22411.09800-red.svg?style=flat)](https://arxiv.org/abs/2411.09800)

# A Probabilistic Framework for Estimating the Modal Age at Death

Manuscript and code for *A Probabilistic Framework for Estimating the Modal Age at Death* by Silvio C. Patricio and Paola Vazquez-Castillo (Interdisciplinary Center on Population Dynamics, University of Southern Denmark).

The method gives, for each age, the posterior probability that it is the modal age at death. Death counts are Poisson given the age-specific death rates and exposures; the posterior of the rates is passed through the life table, and the age with the most life-table deaths is recorded in each posterior draw. Two priors are implemented:

- **independent gamma priors** (the "raw" posterior), with a closed-form posterior of the rates (`Code/mode_posterior.R`);
- a **second-order random-walk (RW2) prior** on the log rates, which lets neighbouring ages share information, computed with a Laplace approximation that integrates over the smoothing strength (`Code/mode_posterior_rw.R`, with the numerical core in C++, `Code/rw_core.cpp`).

The code also contains the three estimators the method is compared with: Gompertz, P-spline and the discretized derivative tests (DDT).

## Structure

```
├── Manuscript/
│   ├── main.tex                  LaTeX source of the paper
│   ├── main.pdf                  compiled paper
│   ├── references.bib
│   ├── figures/                  figures of the paper (written by the code)
│   └── tables/                   LaTeX tables of the paper (written by the code)
│
└── Code/
    ├── 01_simulation.R           Monte Carlo experiment (Section 3)
    ├── 01b_simulation_bootstrap.R  bootstrap intervals of the benchmarks (Section 3)
    ├── 02_figures_simulation.R   figures of Section 3 and Appendix E
    ├── 03_tables_simulation.R    tables of Appendix E
    ├── 04_empirical.R            application to HMD data, women (Section 4)
    ├── 05_figures_empirical.R    figures of Section 4 and Appendix F, Table 1
    ├── 06_laplace_vs_mcmc.R      Laplace approximation vs exact sampling (Appendix D)
    ├── 07_text_numbers.R         every number quoted in the text
    │
    ├── mode_posterior.R          posterior of the mode, independent gamma priors
    ├── mode_posterior_rw.R       posterior of the mode, RW2 prior
    ├── rw_core.cpp               C++ (Rcpp) core of the RW2 posterior and the exact sampler
    ├── mode_gompertz.R           Gompertz estimator
    ├── mode_pspline.R            P-spline estimator
    ├── mode_ddt.R                DDT estimator (Vazquez-Castillo et al. 2024)
    ├── simulation_truth.R        true distributions of the simulation study
    ├── plot_utils.R              plotting helpers for the simulation figures
    ├── utils.R                   data input, life table, parallel computation, saving
    ├── Code.Rproj                RStudio project (sets the working directory)
    │
    ├── data/                     HMD country folders (see below)
    ├── figures/                  figures produced by the scripts
    └── results/                  results and LaTeX tables produced by the scripts
```

## Reproducing the results

Open `Code/Code.Rproj` in RStudio, or start R with `Code/` as working directory, and run the numbered scripts in order:

```r
source("01_simulation.R")          # about 5 minutes on 8 cores
source("01b_simulation_bootstrap.R")  # about 30 minutes on 8 cores
source("02_figures_simulation.R")
source("03_tables_simulation.R")
source("04_empirical.R")           # a few minutes on 8 cores
source("05_figures_empirical.R")
source("06_laplace_vs_mcmc.R")     # about 1 minute
source("07_text_numbers.R")
```

Script 01b needs the output of 01, and scripts 02 and 03 need the output of 01 and 01b; scripts 05 and 07 need the output of 04 (and 07 also that of 01). Script 06 only needs the data. The `Code/results/` and `Code/figures/` folders already contain the outputs used in the paper, so any script can be run on its own.

Each figure and LaTeX table used in the paper is saved both in `Code/figures/` or `Code/results/` and in `Manuscript/figures/` or `Manuscript/tables/`, from where `main.tex` includes it. Every number quoted in the text is recomputed by `07_text_numbers.R` (written to `Code/results/text_numbers.txt`); those of Appendix D are printed by `06_laplace_vs_mcmc.R`.

## Building the manuscript

From `Manuscript/`:

```
pdflatex main
biber main
pdflatex main
pdflatex main
```

or `latexmk -pdf main.tex`. The bibliography uses `biblatex-chicago` with `biber`.

## Requirements

- R (tested with R 4.3).
- R packages: `Rcpp`, `ggplot2`, `patchwork`. The packages `parallel`, `splines` and `stats` come with R.
- A C++ compiler for `Rcpp` (on macOS, the Xcode command line tools; on Windows, Rtools). `rw_core.cpp` is compiled automatically when `mode_posterior_rw.R` is sourced. Without a compiler, the RW2 posterior falls back to a slower R implementation, and `06_laplace_vs_mcmc.R` cannot be run.
- The scripts use forked parallelism (`parallel::mclapply`) on macOS and Linux, with all cores but one. On Windows they run on one core.

## Data

`Code/data/` contains the Human Mortality Database folders of four countries: `DNK` (Denmark), `FRATNP` (France, total population), `JPN` (Japan) and `SWE` (Sweden). The code reads `<code>/STATS/Deaths_1x1.txt` and `<code>/STATS/Exposures_1x1.txt` (single years of age and calendar years) and uses the data for women. The HMD publishes its data under the CC BY 4.0 license. Please cite the database when using them:

> HMD. Human Mortality Database. Max Planck Institute for Demographic Research (Germany), University of California, Berkeley (USA), and French Institute for Demographic Studies (France). Available at www.mortality.org.

The first line of each file gives the date of its HMD revision. Newer versions of the country folders can be dropped into `Code/data/` without changes to the code; the results will then differ slightly from those in the paper. Other countries can be added to `hmd_countries` in `Code/utils.R`.

## License

The code in `Code/` is released under the MIT License (see `LICENSE`). The manuscript in `Manuscript/` is copyright of the authors. The HMD data keep their own license (CC BY 4.0).
