[![arXiv shield](https://img.shields.io/badge/arXiv-2411.09800-red.svg?style=flat)](https://arxiv.org/abs/2411.09800)

# A Probabilistic Framework for Estimating the Modal Age at Death

Manuscript and code for *A Probabilistic Framework for Estimating the Modal Age at Death* by Silvio C. Patricio and Paola Vazquez-Castillo (Interdisciplinary Center on Population Dynamics, University of Southern Denmark).

The method gives, for each age, the posterior probability that it is the modal age at death. Death counts are Poisson given the age-specific death rates and exposures; the posterior of the rates is passed through the life table, and the age with the most life-table deaths is recorded in each posterior draw. Two priors are implemented:

- **independent gamma priors** (the "raw" posterior), with a closed-form posterior of the rates (`Code/mode_posterior.R`);
- a **second-order random-walk (RW2) prior** on the log rates, which lets neighboring ages share information, computed with a Laplace approximation that integrates over the smoothing strength (`Code/mode_posterior_rw.R`, with the numerical core in C++, `Code/rw_core.cpp`).

The code also contains the three estimators the method is compared with: Gompertz, P-spline and the discretized derivative tests (DDT).

The arXiv preprint (v1, 2024) used an earlier version of the method, based on a multinomial model and a Gaussian approximation. That version's code is kept in `legacy/` for reference; everything else in this repository corresponds to the current manuscript in `Manuscript/`.

## Structure

```
├── README.md
├── LICENSE                       MIT license (code)
├── CITATION.cff                  how to cite this repository
├── legacy/
│   └── estimation_algorithm_arxiv_v1.R   code of the arXiv v1 method (multinomial, Gaussian)
│
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
    ├── 06_laplace_vs_mcmc.R      Laplace approximation vs posterior sampling (Appendix D)
    ├── 07_text_numbers.R         every number quoted in the text
    ├── 08_review_checks.R        change 2019-2020 (Table 2, Figure 6) and checks of Appendix G
    │
    ├── mode_posterior.R          posterior of the mode, independent gamma priors
    ├── mode_posterior_rw.R       posterior of the mode, RW2 prior
    ├── rw_core.cpp               C++ (Rcpp) core of the RW2 posterior and the MCMC check
    ├── mode_gompertz.R           Gompertz estimator
    ├── mode_pspline.R            P-spline estimator
    ├── mode_ddt.R                DDT estimator (Vazquez-Castillo et al. 2024)
    ├── simulation_truth.R        true distributions of the simulation study
    ├── plot_utils.R              plotting helpers for the simulation figures
    ├── utils.R                   data input, life table, parallel computation, saving
    ├── Code.Rproj                RStudio project (sets the working directory)
    │
    ├── data/                     HMD deaths and exposures, four countries (see below)
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
source("06_laplace_vs_mcmc.R")     # a few minutes
source("07_text_numbers.R")
source("08_review_checks.R")       # about 30 minutes on 8 cores
```

Script 01b needs the output of 01, and scripts 02 and 03 need the output of 01 and 01b; scripts 05 and 07 need the output of 04 (and 07 also that of 01); script 08 needs the output of 01, 01b and 04. Script 06 only needs the data. Part 4 of script 08 (about 30 minutes) can be skipped with `Sys.setenv(MODE_SKIP_PSPLINE = "1")`. The `Code/results/` and `Code/figures/` folders already contain the outputs used in the paper, so any script can be run on its own.

Each figure and LaTeX table used in the paper is saved both in `Code/figures/` or `Code/results/` and in `Manuscript/figures/` or `Manuscript/tables/`, from where `main.tex` includes it. Every number quoted in the text is recomputed by `07_text_numbers.R` (written to `Code/results/text_numbers.txt`); those of Appendix D are printed by `06_laplace_vs_mcmc.R`, and those of Table 2 and Appendix G are written by `08_review_checks.R` to `Code/results/review_checks/review_numbers.txt`.

## Building the manuscript

From `Manuscript/`:

```
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

or `latexmk -pdf main.tex`. The bibliography uses `natbib` with the `apalike` style. With TinyTeX in R, `tinytex::latexmk("main.tex")` also installs any missing LaTeX packages.

## Requirements

- R (tested with R 4.3).
- R packages: `Rcpp`, `ggplot2`, `patchwork`. The packages `parallel`, `splines` and `stats` come with R.
- A C++ compiler for `Rcpp` (on macOS, the Xcode command line tools; on Windows, Rtools). `rw_core.cpp` is compiled automatically when `mode_posterior_rw.R` is sourced. Without a compiler, the RW2 posterior falls back to a slower R implementation, and `06_laplace_vs_mcmc.R` cannot be run.
- The scripts use forked parallelism (`parallel::mclapply`) on macOS and Linux, with all cores but one. On Windows they run on one core.

## Data

`Code/data/` contains the Human Mortality Database files used in the paper for four countries: `DNK` (Denmark), `FRATNP` (France, total population), `JPN` (Japan) and `SWE` (Sweden). For each country, only `<code>/STATS/Deaths_1x1.txt` and `<code>/STATS/Exposures_1x1.txt` (single years of age and calendar years) are included; the code uses the data for women. The data were downloaded on October 28, 2026; the first line of each file gives the date of its HMD revision. The HMD publishes its data under the CC BY 4.0 license. Please cite the database when using them:

> HMD. Human Mortality Database. Max Planck Institute for Demographic Research (Germany), University of California, Berkeley (USA), and French Institute for Demographic Studies (France). Available at www.mortality.org.

Newer versions of these files, or the full HMD country folders, can be dropped into `Code/data/` without changes to the code; the results will then differ slightly from those in the paper. Other countries can be added to `hmd_countries` in `Code/utils.R`.

## Citation

If you use this code, please cite the paper (see `CITATION.cff`):

> Patricio, S. C. and Vazquez-Castillo, P. A Probabilistic Framework for Estimating the Modal Age at Death. arXiv:2411.09800.

## License

The code in `Code/` is released under the MIT License (see `LICENSE`). The manuscript in `Manuscript/` is copyright of the authors. The HMD data keep their own license (CC BY 4.0).
