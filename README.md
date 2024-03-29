# ChenPhelanMPStability

Replication code for
"Should Monetary Policy Target Financial Stability?" by
William Chen and Gregory Phelan

The code in this repository produces all figures and tables found
in the paper, plus some additional unreported results.
The code was run on MacOS; it has not been tested
on other operating systems.
Some results require a long computation time, e.g. the
optimal monetary policy rule. The scripts for these latter
results are also provided, but we recommend running
those scripts on a server in the background.

If the user is interested in using this repository's code, then
looking at `fedput_calibration.m` will be a good starting point.
The key function for solving equilibrium is `eqm_functions/geteqm.m`.
That function's source code will be helpful in learning how the
rest of the source code works since `geteqm.m` is intended to be
the user-level interface with the underlying least-squares problem.
Please file an issue if you have any questions.

## Directory Structure

* MATLAB scripts in top-level of directory: produce results for paper's main body

* `eqm_functions/`: functions used to compute equilibrium and welfare

* `parameters_specs/`: scripts that initialize economic and numerical parameters.
The code assumes that the default interpreter for plot text is LaTeX,
but the parameter specfiles automatically sets this up.

* `save/`: input data for data moments and saved output

* `figures/`: folder for figures

* `appendix/`: scripts for calculating figures and tables in the Appendix

* `batch_jobs/`: scripts for running batch jobs

## Replicating Tables and Figures in Paper's Main Body

Below, we report which script generates which
figure/table in the main body of the paper. This section assumes output has already been generated
by the batch scripts, and the user only wants to produce tables and figures.
In order to save figures to the `figures/` folder, the variable
`save_plots` near the top of every script must be set to true (or 1).
If none of the scripts have been edited, then
running `run_all.m` with `main_paper = 1`
will run all of the scripts below.

- Table 2: `fedput_calibration.m`
- Figure 1: `fedput_calibration.m`
- Figure 2: `optnomap_constpi_friedman.m`
- Table 3: `optnomap_comparison.m`
- Figure 3: `optnomap_constpi_friedman.m`
- Figure 4: `lc_icc_example_plots.m`
- Figure 5: `OptMP_vary_alpha.m`

## Replicating Tables and Figures in the Appendix

Below, we report which script generates which
figure/table in Appendices C and E.
The "current working directory" for these scripts
is `appendix/`. This section assumes output has already been generated
by the batch scripts, and the user only wants to produce tables and figures.
In order to save figures to the `figures/` folder, the variable
`save_plots` near the top of every script must be set to true (or 1).
We do not provide the code for Appendix D.
If none of the scripts have been edited, then
running `run_all.m` with `appendix = 1`
will run all of the scripts below.

- Table 4: `fedput_calibration_constpi.m`
- Figure 6: `optnomap_vary_lambda.m`
- Table 3 welfare comparisons for different lambda: `optnomap_comparison_lambda.m`.
- Figure 7: `optnomap_vary_iput.m`
- Figure 8: `optnomap_vary_ilaw.m`
- Figure 9: `optnomap_vary_strike.m`
- Figure 10: `optnomap_constrates.m`
- Figure 11: `../optnomap_constpi_friedman.m`
- Table 5: `map_comparison_well_targeted.m`
- Figure 12: `stability_fraction_03.m` and `stability_fraction_07.m`
- Figure 13: `lc_poorly_targeted_example_plots.m`
- Table 6: `map_comparison_poorly_targeted.m`
- Figure 15: `bank_misalloc/fedput_calibration.m`
- Figure 16: `bank_misalloc/optnomap_vs_no_bank_misalloc.m`

## Computing Optimal MP Rules

The following steps will calculate the various optimal
MP rules in the paper using batch jobs.

1. Run `max_expV_fedput_4args.m` to compute OptNoMaP.
2. Run all of the `max_expV_fedput_4args_icc*.m` scripts to calculate OptMP as a function of alpha
   for `OptMP_vary_alpha.m`
3. Run both of the `max_expV_fedput_4args_lambda*.m` scripts to calculate OptNoMaP
   for different lambda values.
4. Run `bank_misalloc_max_expV_fedput_4args.m` to calculate OptNoMaP
   when banks misallocate capital.

## Miscellaneous Scripts

- We provide example scripts that
  solve the model with exogenous equity injections in `appendix/equity_injections/`

- The script `appendix/calculate_avg_gdp_growth_vol.m` calculates the average volatility of annual GDP growth
  using data from the Penn World Tables.
