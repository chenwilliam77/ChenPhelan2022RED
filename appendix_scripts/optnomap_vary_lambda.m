% Create plots of the optimal monetary policy rule with no MaP
% and the associated value function
% under different values of lambda
%
% Written by William Chen and Gregory Phelan, Feb. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

parameters_file = 'calibration_parameters';
figurespath = [pwd(), '/../figures/optnomap_vary_lambda/'];
save_plots = 1; % Set to 1 if you want to save plots
find_constpi_rule = 1; % Set to 1 to calculate the constant inflation MP rule, otherwise load from memory
constpi_path = [pwd(), '/../save/output_data/fedput_calibration_constpi/']; % save path for constant inflation MP rule
load(['../save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat']);

leg_text = {'OptNoMaP', 'Inflation Targeting', '1/2 $\lambda$', '$\lambda = 0$'};

%% Compute equilibrium
% Fed put
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
fp_s = s; % copy the struct
[fp_solution, fp_summary_stats] = geteqm(fp_s);

% constant inflation
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.color = 'b';
s.nomrule = 2; % turn on interpolant for MP rule
if find_constpi_rule
    [constpi_s, converged] = ...
        find_constant_inflation_mp_rule(s, 'damping', 0.7, 'tol', 5e-3);
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
[constpi_solution, constpi_summary_stats] = geteqm(constpi_s);

% lambda = 1/2 calibration
eval(parameters_file);
load(['../save/output_data/max_expV_fedput_4args_lambda_half/fmincon_run_interior-point_08-Mar-2022.mat']);
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.color = [0.4940, 0.1840, 0.5560];
half_s = s;
[half_solution, half_summary_stats] = geteqm(half_s);

% lambda = 0
eval(parameters_file);
load(['../save/output_data/max_expV_fedput_4args_lambda0/fmincon_run_interior-point_08-Mar-2022.mat']);
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.color = [0.3010, 0.7450, 0.9330];
zero_s = s;
[zero_solution, zero_summary_stats] = geteqm(zero_s);

% Make initial plots
plot_results(fp_solution, fp_s);
plot_results(constpi_solution, constpi_s);
plot_results(half_solution, half_s);
plot_results(zero_solution, zero_s);

%% Edits to plots
fignum = [32, 56];
save_fn = {'nominal_interest', 'welfare_cons'};

% Add legends and save if needed
for idx = 1:length(fignum)
    figure(fignum(idx));
    if fignum(idx) == 32
        legend(leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');
    end
end

% Save plots
if save_plots
    if ~isfolder(figurespath)
        mkdir(figurespath);
    end

    for i = 1:length(fignum)
        f = figure(fignum(i));
        saveas(f, [figurespath, '/', save_fn{i}, '.', s.graphics_format]);
    end
end
