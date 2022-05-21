# ChenPhelanMPStability

Replication code for
"Should Monetary Policy Target Financial Stability?" by
William Chen and Gregory Phelan

The code in this repository produces the main figures and tables.
Some results require a long computation time, e.g. the
optimal monetary policy rule. The scripts for these latter
results are also provided, but we recommend running
those scripts on a server in the background.

If the user is interested in using this repository's code, then
looking at `fedput_only.m` will be a good starting point.
The key function for solving equilibrium is `eqm_functions/geteqm.m`.
That function's source code will be helpful in learning how the
rest of the source code works since `geteqm.m` is intended to be
the user-level interface with the underlying least-squares problem.
Please file an issue if you have any questions.

## Directory Structure

* MATLAB scripts in top-level of directory: produce results for paper's main body

* `eqm_functions/`: functions used to compute equilibrium and welfare for the
baseline model and the extension to CES production

* `parameter_specs/`: scripts that initialize economic and numerical parameters.
MAKE A NOTE HERE THAT THE CODE ASSUMES THE DEFAULT INTERPRETER FOR PLOT TEXT IS LATEX.

* `save/`: input data for data moments and saved output

* `figures/`: folder for figures

* `appendix/`: scripts for calculating figures and tables in the Appendix

## Tables and Figures Replication Instructions

The following instructions assumes output has already been generated
and the user only wants to produce tables and figures.

1. Run `make_paper_figs.m` in the top-level, `1good/`, and `linear_invst/`;
   and `dominant_strat.m` in `1good/` and `linear_invst/` to produce all the
   main tables and figures reported in the main text. Note that the numbers in
   tables were manually inputted, so these scripts will only generate the
   numbers reported, not the actual tables themselves.

2. The script `make_paper_figs` will sketch out the welfare surface
   for the AE and EME as a function of their leverage constraint policies.
   These welfare surfaces visually show the numerical result that the
   zero constraint is a dominant strategy.
   To further check the Nash result, run `analyze_check_zero.m` and
   `analyze_iternash.m`, which analyze output from the two other
   approaches for checking the Nash result.

3. Run `compare_stat_dist.m` to generate the counterfactual stationary distributions
   in Appendix B.3.

4. Run `bargaining_example.m` to obtain the numbers reported in the bargaining
   example of Appendix D.1.

5. Run `analyze_proddiff_robustness.m` and `analyze_risk_aversion_robustness.m`
   to obtain the numbers reported in Appendices E.1 - E.3.

6. Run `make_paper_ces_tbls_figs.m` to obtain the numbers reported in Appendix E.6
   on the extension to CES production.

7. Run `compare_ce_nash_specific_init_conds.m` to obtain the comparison
   of welfare in the competitive equilibrium to Nash at specific initial conditions.
   These numbers are discussed

## Output Replication Instructions

To re-generate all the saved output from scratch, follow these instructions.

1. Run `numerical_example.m` with `use_saved_guess = 0` to generate the
   converged solution in `data/gamA1p05_gamB2.mat`.

2. Run `compute_eqm_for_benchmarks_etc.m` and `compute_ces_eqms.m` to generate
   some output that is needed by `1good/` and `linear_invst/`

3. Change directory to `hpcc/`.

4. Run `hpcc/plot_best_response.m` to sketch out the best responses
   and provide a coarse grid of initial guesses for the other scripts.

5. Run `plot_stat_pdf.m`, `get_best_response.m`, `sparse_get_br.m`, `check_zero.m`, and `iter_nash.m`
   in `hpcc/`.

6. Run `analyze_policy.m` twice after `get_best_response.m` and `sparse_get_br.m` finish, once each
   for the output of the latter two scripts.

7. To check the robustness results for different parameters, re-run steps 4-6 but either adding lines
   to the scripts to change parameters or creating a new parameter script instead of `parameters/baseline_parameters.m`.
