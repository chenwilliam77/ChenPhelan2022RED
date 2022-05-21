% Welfare and stability numbers comparing monetary policy rules
% to OptNoMaP when lambda is half its calibrated value
%
% Written by William Chen, Sep. 2020

close all;
clear;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

parameters_file = 'calibration_parameters';
load(['../save/output_data/max_expV_fedput_4args_lambda_half/fmincon_run_interior-point_08-Mar-2022.mat']);

%% Rule 1 - Optimal Fed Put
n_comps = 5;
rule1_stats.welfare = zeros(n_comps, 1);
rule1_stats.stab_frac = zeros(n_comps, 1);
rule1_stats.stab_frac_pct50 = zeros(n_comps, 1);
rule1_sols = cell(n_comps, 1);
rule1_sumstats = cell(n_comps, 1);

eval(parameters_file);
s.lambda = s.lambda / 2;
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;

% Best rule
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);

L1s = s;
[rule1_sols{1}, rule1_sumstats{1}] = geteqm(L1s);

% OptNoMaP
rule1_stats.welfare(1) = 0;
rule1_stats.stab_frac(1) = 0;

% i = iput constant
L1s.ilaw = best_rule(1);
[rule1_sols{2}, rule1_sumstats{2}] = geteqm(L1s);
rule1_stats.welfare(2) = consumption_equiv(rule1_sumstats{2}.ergodic.ExpV, rule1_sumstats{1}.ergodic.ExpV, L1s.r);
rule1_stats.stab_frac(2) = rule1_sumstats{2}.ergodic.stab_frac - rule1_sumstats{1}.ergodic.stab_frac;
rule1_stats.stab_frac_pct50(2) = rule1_sumstats{2}.ergodic.stab_frac_pct50 - rule1_sumstats{1}.ergodic.stab_frac_pct50;

% i = ilaw constant
L1s.ilaw = best_rule(2);
L1s.iput = best_rule(2);
[rule1_sols{3}, rule1_sumstats{3}] = geteqm(L1s);
rule1_stats.welfare(3) = consumption_equiv(rule1_sumstats{3}.ergodic.ExpV, rule1_sumstats{1}.ergodic.ExpV, L1s.r);
rule1_stats.stab_frac(3) = rule1_sumstats{3}.ergodic.stab_frac - rule1_sumstats{1}.ergodic.stab_frac;
rule1_stats.stab_frac_pct50(3) = rule1_sumstats{3}.ergodic.stab_frac_pct50 - rule1_sumstats{1}.ergodic.stab_frac_pct50;

% i = Friedman with lambda > 0
friedman_rate_lambda_pos = compute_first_best_rate(L1s);
L1s.ilaw = friedman_rate_lambda_pos;
L1s.iput = friedman_rate_lambda_pos;
[rule1_sols{4}, rule1_sumstats{4}] = geteqm(L1s);
rule1_stats.welfare(4) = consumption_equiv(rule1_sumstats{4}.ergodic.ExpV, rule1_sumstats{1}.ergodic.ExpV, L1s.r);
rule1_stats.stab_frac(4) = rule1_sumstats{4}.ergodic.stab_frac - rule1_sumstats{1}.ergodic.stab_frac;
rule1_stats.stab_frac_pct50(4) = rule1_sumstats{4}.ergodic.stab_frac_pct50 - rule1_sumstats{1}.ergodic.stab_frac_pct50;

% Inflation Targeting
[constpi_s, converged] = ...
    find_constant_inflation_mp_rule(L1s, 'damping', 0.7, 'tol', 5e-3);
if ~converged
    disp('Inflation targeting rule not found! Figure out what is wrong');
    keyboard;
end
L1s = constpi_s;
[rule1_sols{5}, rule1_sumstats{5}] = geteqm(L1s);
rule1_stats.welfare(5) = consumption_equiv(rule1_sumstats{5}.ergodic.ExpV, rule1_sumstats{1}.ergodic.ExpV, L1s.r);
rule1_stats.stab_frac(5) = rule1_sumstats{5}.ergodic.stab_frac - rule1_sumstats{1}.ergodic.stab_frac;
rule1_stats.stab_frac_pct50(5) = rule1_sumstats{5}.ergodic.stab_frac_pct50 - rule1_sumstats{1}.ergodic.stab_frac_pct50;

%% Show results!
fprintf('(ExpV, Pr. Distress, Pr. Crisis) for Rule 1:     (%.2f, %.2f%%, %.2f%%)\n', ...
        [rule1_sumstats{1}.ergodic.ExpV, 100 - rule1_sumstats{1}.ergodic.stab_frac, ...
        100 - rule1_sumstats{1}.ergodic.stab_frac_pct50]);

fprintf('\n')
fprintf('\n')
fprintf('\n')

policy_names = {'OptNoMaP', 'i = optimal iput', 'i = optimal ilaw', 'i = Friedman with pos. lambda', 'inflation targeting'}; % 'higher LAW', 'lower LAW', 'higher Put'};
for i = 1:n_comps
    fprintf(['(ExpV, Pr. Distress, Pr. Crisis) Differences for Rule 1, ', policy_names{i}, ':     (%.2f%%, %.2f%%, %.2f%%)\n'], ...
            [100 * rule1_stats.welfare(i), -rule1_stats.stab_frac(i), ...
            -rule1_stats.stab_frac_pct50(i)]);
end

fprintf('\n')
