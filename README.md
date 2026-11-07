This repository contains data and analysis scripts accompanying the paper
[*Modeling Emerging Online Language: A Framework for Analyzing Neologism Diffusion Over Time*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7532402),
currently available as an SSRN preprint.

This folder contains annual Reddit and YouTube time series and a small R workflow
to fit Bass diffusion curves, export parameter estimates, and produce descriptive
summaries and figures.

## Requirements

Use R 4.2 or later and install the fitting and plotting packages:

```r
install.packages(c("minpack.lm", "ggplot2", "ggrepel", "patchwork", "paletteer", "scales"))
```

Tested with R 4.6.1. Package versions are recorded in `results/summary.txt`.

## Run everything

From a terminal in this folder:

```sh
Rscript run_all.R
```

Alternatively, set the working directory to this folder and run:

```r
source("run_all.R")
```

Results are written to `results/`.

## Included data

The two CSVs in `data/` contain `year`, `neologism`, and `count`.
The platform is identified by the filename. Counts are annual term-comment
matches for 2006--2024 obtained from the original dataset.

`neologism_categories.csv` provides the semantic categories and display labels used in the plots. 
The `label_to_show` column contains labels only for non-offensive neologisms. 
Offensive terms are left unlabelled in the figures.
For further information about the extraction pipeline, you are referred to the 
methods section of the [paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7532402).

For further information about the extraction pipeline, you are referred to the 
methods section of the [paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7532402).

## Fitting procedure

The default analysis retains trajectories with at least 500 total uses and six
annual observations. Observation duration includes intervening zero-count years.
These settings and the study period can be changed in
the `settings` list near the beginning of `run_all.R`.

The model is fitted to cumulative counts:

```text
N(t) = m * (1 - exp(-(p + q) * t)) / (1 + (q / p) * exp(-(p + q) * t))
```

The first annual cumulative observation is at `t = 1`, with `N(0) = 0`.
The parameters `m`, `p`, and `q` are estimated jointly and the lower bound on `m`
is the largest observed cumulative count.
The original OLS initialization and nonlinear fitting/retry procedure are
implemented in `bass_model.R`.

Approximate parameter intervals use the fitted covariance matrix on the log scale.

## Outputs

- `trajectory_summary.csv`: inclusion decisions and fitting status for every input trajectory.
- `bass_parameters.csv`: fitted parameters, approximate intervals, and fitting errors.
- `platform_summary.csv`: numbers of trajectories and median parameter/error summaries.
- `summary.txt`: settings, summary tables, figure descriptions, and software versions.
- `figures/example_fits.png`: fitting of four example curves.
- `figures/bass_parameters.png`: category-coloured parameter plots and cross-platform comparisons, as in the paper. The parameter-plane panel includes fits with R-squared at least 0.95.
- `figures/fit_error.png`: NRMSE distributions by platform.
