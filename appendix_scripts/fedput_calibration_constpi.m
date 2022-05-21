% This script is primarily used to report moments under
% Calibration and inflation targeting.
%
% Written by William Chen and Greg Phelan, Mar. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

parameters_file = 'calibration_parameters';

figurespath = [pwd(), '/../figures/fedput_calibration_constpi_friedman/'];
save_plots = 1; % Set to 1 if you want to save plots
find_constpi_rule = 1; % Set to 1 to calculate the constant inflation MP rule, otherwise load from memory
constpi_path = [pwd(), '/../save/output_data/fedput_calibration_constpi/']; % save path for constant inflation MP rule
leg_text = {'Calibration', 'Inflation Targeting', 'Friedman'};

%% Compute equilibrium
% Fed put
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.crisis_output_losses = 1;
fp_s = s; % copy the struct
rng(1793);
[fp_solution, fp_summary_stats] = geteqm(fp_s);
fp_s.crisis_output_losses = 0;
s.crisis_output_losses = 0;

% constant inflation
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.color = 'b';
s.nomrule = 2; % turn on interpolant for MP rule
s.verbose = 0;
if find_constpi_rule
    [constpi_s, converged] = find_constant_inflation_mp_rule(s, 'tol', 5e-3, 'damping', 0.7);
    if converged
        constpi_nomrule_mp_interp = constpi_s.nomrule_mp_interp;
        save([constpi_path, 'fedput_calibration_constpi_nomrule_mp_interp.mat'], 'constpi_nomrule_mp_interp');
    else
        error('Convergence not achieved for finding the constant inflation rule');
    end
else
    constpi_s = s;
    load([constpi_path, 'fedput_calibration_constpi_nomrule_mp_interp.mat'], 'constpi_nomrule_mp_interp');
    constpi_s.nomrule_mp_interp = constpi_nomrule_mp_interp;
end
constpi_s.crisis_output_losses = 1;
constpi_s.verbose = 1;
rng(1793);
[constpi_solution, constpi_summary_stats] = geteqm(constpi_s);
constpi_s.crisis_output_losses = 0;
