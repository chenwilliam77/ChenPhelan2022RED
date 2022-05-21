% Welfare and stability numbers for well-targeted macroprudential policy (ICC)
% under different monetary policy rules
%
% Written by William Chen and Gregory Phelan, Mar. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

parameters_file = 'lc_icc_calibration_parameters';

% Settings for what to do
N_mp_rules = 5;
alpha_vec = [Inf; 8; 6];
day_cells = {''; '08'; '08'};

%% Rule 1
% Load Optimal Fed Put
load('../save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat');
disp('Running rule 1');

% Set storage containers
rule1_stats.welfare = zeros(N_mp_rules, 1);
rule1_stats.stab_frac = zeros(N_mp_rules, 1);
rule1_sols = cell(N_mp_rules, 1);
rule1_sumstats = cell(N_mp_rules, 1);

eval(parameters_file);
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
s.leverage_constraint.type = 'icc';
s.leverage_constraint.alpha = alpha_vec(1);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
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
[constpi_s, converged] = find_constant_inflation_mp_rule(L1s, 'tol', 5e-3, 'damping', 0.7);
if ~converged
    disp('Inflation targeting rule not found! Figure out what is wrong');
    keyboard;
end
L1s = constpi_s;
[rule1_sols{5}, rule1_sumstats{5}] = geteqm(L1s);
rule1_stats.welfare(5) = consumption_equiv(rule1_sumstats{5}.ergodic.ExpV, rule1_sumstats{1}.ergodic.ExpV, L1s.r);
rule1_stats.stab_frac(5) = rule1_sumstats{5}.ergodic.stab_frac - rule1_sumstats{1}.ergodic.stab_frac;
rule1_stats.stab_frac_pct50(5) = rule1_sumstats{5}.ergodic.stab_frac_pct50 - rule1_sumstats{1}.ergodic.stab_frac_pct50;

%% Rule 2
disp('Running rule 2');

rule2_stats.welfare = zeros(N_mp_rules, 1);
rule2_stats.stab_frac = zeros(N_mp_rules, 1);
rule2_sols = cell(N_mp_rules, 1);
rule2_sumstats = cell(N_mp_rules, 1);
load(['../save/output_data/max_expV_fedput_4args_icc' num2str(alpha_vec(2)) ...
    '/fmincon_run_interior-point_' day_cells{2} '-Mar-2022.mat']);

eval(parameters_file);
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
s.leverage_constraint.type = 'icc';
s.leverage_constraint.alpha = alpha_vec(2);
s.plot_results = 0;
s.use_title = 0;
s.color = 'b';
s.verbose = 0;
L2s = s;
[rule2_sols{1}, rule2_sumstats{1}] = geteqm(L2s);

% OptMP w/ icc
rule2_stats.welfare(1) = consumption_equiv(rule2_sumstats{1}.ergodic.ExpV, rule2_sumstats{1}.ergodic.ExpV, L2s.r);
rule2_stats.stab_frac(1) = rule2_sumstats{1}.ergodic.stab_frac - rule2_sumstats{1}.ergodic.stab_frac;

% i = iput constant
L2s.ilaw = best_rule(1);
[rule2_sols{2}, rule2_sumstats{2}] = geteqm(L2s);
rule2_stats.welfare(2) = consumption_equiv(rule2_sumstats{2}.ergodic.ExpV, rule2_sumstats{1}.ergodic.ExpV, L2s.r);
rule2_stats.stab_frac(2) = rule2_sumstats{2}.ergodic.stab_frac - rule2_sumstats{1}.ergodic.stab_frac;
rule2_stats.stab_frac_pct50(2) = rule2_sumstats{2}.ergodic.stab_frac_pct50 - rule2_sumstats{1}.ergodic.stab_frac_pct50;

% i = ilaw constant
L2s.ilaw = best_rule(2);
L2s.iput = best_rule(2);
[rule2_sols{3}, rule2_sumstats{3}] = geteqm(L2s);
rule2_stats.welfare(3) = consumption_equiv(rule2_sumstats{3}.ergodic.ExpV, rule2_sumstats{1}.ergodic.ExpV, L2s.r);
rule2_stats.stab_frac(3) = rule2_sumstats{3}.ergodic.stab_frac - rule2_sumstats{1}.ergodic.stab_frac;
rule2_stats.stab_frac_pct50(3) = rule2_sumstats{3}.ergodic.stab_frac_pct50 - rule2_sumstats{1}.ergodic.stab_frac_pct50;

% i = Friedman with lambda > 0
friedman_rate_lambda_pos = compute_first_best_rate(L2s);
L2s.ilaw = friedman_rate_lambda_pos;
L2s.iput = friedman_rate_lambda_pos;
[rule2_sols{4}, rule2_sumstats{4}] = geteqm(L2s);
rule2_stats.welfare(4) = consumption_equiv(rule2_sumstats{4}.ergodic.ExpV, rule2_sumstats{1}.ergodic.ExpV, L2s.r);
rule2_stats.stab_frac(4) = rule2_sumstats{4}.ergodic.stab_frac - rule2_sumstats{1}.ergodic.stab_frac;
rule2_stats.stab_frac_pct50(4) = rule2_sumstats{4}.ergodic.stab_frac_pct50 - rule2_sumstats{1}.ergodic.stab_frac_pct50;

% Inflation Targeting
[constpi_s, converged] = find_constant_inflation_mp_rule(L2s, 'damping', 0.7, 'thinning', 5);
if ~converged
    disp('Inflation targeting rule not found! Figure out what is wrong');
    keyboard;
end
L2s = constpi_s;
[rule2_sols{5}, rule2_sumstats{5}] = geteqm(L2s);
rule2_stats.welfare(5) = consumption_equiv(rule2_sumstats{5}.ergodic.ExpV, rule2_sumstats{1}.ergodic.ExpV, L2s.r);
rule2_stats.stab_frac(5) = rule2_sumstats{5}.ergodic.stab_frac - rule2_sumstats{1}.ergodic.stab_frac;
rule2_stats.stab_frac_pct50(5) = rule2_sumstats{5}.ergodic.stab_frac_pct50 - rule2_sumstats{1}.ergodic.stab_frac_pct50;

%% Rule 3
disp('Running rule 3');

rule3_stats.welfare = zeros(N_mp_rules, 1);
rule3_stats.stab_frac = zeros(N_mp_rules, 1);
rule3_sols = cell(N_mp_rules, 1);
rule3_sumstats = cell(N_mp_rules, 1);
load(['../save/output_data/max_expV_fedput_4args_icc' num2str(alpha_vec(3)) ...
    '/fmincon_run_interior-point_' day_cells{3} '-Mar-2022.mat']);

eval(parameters_file);
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
s.leverage_constraint.type = 'icc';
s.leverage_constraint.alpha = alpha_vec(3);
s.plot_results = 0;
s.use_title = 0;
s.color = 'r';
s.verbose = 0;
L3s = s;
[rule3_sols{1}, rule3_sumstats{1}] = geteqm(L3s);

% OptMP w/ icc
rule3_stats.welfare(1) = consumption_equiv(rule3_sumstats{1}.ergodic.ExpV, rule3_sumstats{1}.ergodic.ExpV, L3s.r);
rule3_stats.stab_frac(1) = rule3_sumstats{1}.ergodic.stab_frac - rule3_sumstats{1}.ergodic.stab_frac;

% i = iput constant
L3s.ilaw = best_rule(1);
[rule3_sols{2}, rule3_sumstats{2}] = geteqm(L3s);
rule3_stats.welfare(2) = consumption_equiv(rule3_sumstats{2}.ergodic.ExpV, rule3_sumstats{1}.ergodic.ExpV, L3s.r);
rule3_stats.stab_frac(2) = rule3_sumstats{2}.ergodic.stab_frac - rule3_sumstats{1}.ergodic.stab_frac;
rule3_stats.stab_frac_pct50(2) = rule3_sumstats{2}.ergodic.stab_frac_pct50 - rule3_sumstats{1}.ergodic.stab_frac_pct50;

% i = ilaw constant
L3s.ilaw = best_rule(2);
L3s.iput = best_rule(2);
[rule3_sols{3}, rule3_sumstats{3}] = geteqm(L3s);
rule3_stats.welfare(3) = consumption_equiv(rule3_sumstats{3}.ergodic.ExpV, rule3_sumstats{1}.ergodic.ExpV, L3s.r);
rule3_stats.stab_frac(3) = rule3_sumstats{3}.ergodic.stab_frac - rule3_sumstats{1}.ergodic.stab_frac;
rule3_stats.stab_frac_pct50(3) = rule3_sumstats{3}.ergodic.stab_frac_pct50 - rule3_sumstats{1}.ergodic.stab_frac_pct50;

% i = Friedman with lambda > 0
friedman_rate_lambda_pos = compute_first_best_rate(L3s);
L3s.ilaw = friedman_rate_lambda_pos;
L3s.iput = friedman_rate_lambda_pos;
[rule3_sols{4}, rule3_sumstats{4}] = geteqm(L3s);
rule3_stats.welfare(4) = consumption_equiv(rule3_sumstats{4}.ergodic.ExpV, rule3_sumstats{1}.ergodic.ExpV, L3s.r);
rule3_stats.stab_frac(4) = rule3_sumstats{4}.ergodic.stab_frac - rule3_sumstats{1}.ergodic.stab_frac;
rule3_stats.stab_frac_pct50(4) = rule3_sumstats{4}.ergodic.stab_frac_pct50 - rule3_sumstats{1}.ergodic.stab_frac_pct50;

% Inflation Targeting
[constpi_s, converged] = find_constant_inflation_mp_rule(L3s, 'damping', 0.7, 'thinning', 5);
if ~converged
    disp('Inflation targeting rule not found! Figure out what is wrong');
    keyboard;
end
L3s = constpi_s;
[rule3_sols{5}, rule3_sumstats{5}] = geteqm(L3s);
rule3_stats.welfare(5) = consumption_equiv(rule3_sumstats{5}.ergodic.ExpV, rule3_sumstats{1}.ergodic.ExpV, L3s.r);
rule3_stats.stab_frac(5) = rule3_sumstats{5}.ergodic.stab_frac - rule3_sumstats{1}.ergodic.stab_frac;
rule3_stats.stab_frac_pct50(5) = rule3_sumstats{5}.ergodic.stab_frac_pct50 - rule3_sumstats{1}.ergodic.stab_frac_pct50;

%% Show results!
fprintf('(ExpV, Pr. Distress, Pr. Crisis) for Rule 1:     (%.2f, %.2f%%, %.2f%%)\n', ...
        [rule1_sumstats{1}.ergodic.ExpV, 100 - rule1_sumstats{1}.ergodic.stab_frac, ...
        100 - rule1_sumstats{1}.ergodic.stab_frac_pct50]);
fprintf('(ExpV, Pr. Distress, Pr. Crisis) for Rule 2:     (%.2f, %.2f%%, %.2f%%)\n', ...
        [rule2_sumstats{1}.ergodic.ExpV, 100 - rule2_sumstats{1}.ergodic.stab_frac, ...
        100 - rule2_sumstats{1}.ergodic.stab_frac_pct50]);
fprintf('(ExpV, Pr. Distress, Pr. Crisis) for Rule 3:     (%.2f, %.2f%%, %.2f%%)\n', ...
        [rule3_sumstats{1}.ergodic.ExpV, 100 - rule3_sumstats{1}.ergodic.stab_frac, ...
        100 - rule3_sumstats{1}.ergodic.stab_frac_pct50]);

fprintf('\n')
fprintf('\n')
fprintf('\n')

policy_names = {'OptMP', 'i = optimal iput', 'i = optimal ilaw', 'i = Friedman with pos. lambda', 'inflation targeting'}; % 'higher LAW', 'lower LAW', 'higher Put'};
for i = 1:N_mp_rules
    fprintf(['(ExpV, Pr. Distress, Pr. Crisis) Differences for Rule 1, ', policy_names{i}, ':     (%.2f%%, %.2f%%, %.2f%%)\n'], ...
            [100 * rule1_stats.welfare(i), -rule1_stats.stab_frac(i), ...
            -rule1_stats.stab_frac_pct50(i)]);
end

fprintf('\n')
fprintf('\n')
fprintf('\n')

for i = 1:N_mp_rules
    fprintf(['(ExpV, Pr. Distress, Pr. Crisis) Differences for Rule 2, ', policy_names{i}, ':     (%.2f%%, %.2f%%, %.2f%%)\n'], ...
            [100 * rule2_stats.welfare(i), -rule2_stats.stab_frac(i), ...
            -rule2_stats.stab_frac_pct50(i)]);
end

fprintf('\n')
fprintf('\n')
fprintf('\n')

for i = 1:N_mp_rules
    fprintf(['(ExpV, Pr. Distress, Pr. Crisis) Differences for Rule 3, ', policy_names{i}, ':     (%.2f%%, %.2f%%, %.2f%%)\n'], ...
            [100 * rule3_stats.welfare(i), -rule3_stats.stab_frac(i), ...
            -rule3_stats.stab_frac_pct50(i)]);
end

fprintf('\n')
fprintf('\n')
fprintf('\n')
