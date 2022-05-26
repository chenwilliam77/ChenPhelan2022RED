% Create plots of equilibrium dynamics when using the
% Calibration parameters
%
% Written by William Chen and Gregory Phelan, Jan. 2022

clear;
close all;
addpath parameters_specs
addpath eqm_functions

parameters_file = 'calibration_parameters';
save_plots = 1; % Set to 1 if you want to save plots
figurespath = [pwd(), '/figures/fedput_calibration/'];

%% Compute equilibrium
% Fed put
eval(parameters_file);
s.plot_results = 1;
s.use_title = 0;
s.crisis_output_losses = 1; % report crisis output losses in the summary statistics
rng(1793);
[solution, summary_stats] = geteqm(s);
s.crisis_output_losses = 0;

%% Edits to plots
fignum = [1, 2, 22, 23, 6, 7, 8, 9, 10, 11, 12, 16, 25, 26, 32, 33, 35, ...
    40, 41, 42, 43, 56, 58]; % figures to save/edit in loops
save_fn = {'leverage', 'mu_eta', 'sigma_eta', 'psi', 'omega', 'sigma_Ch', 'sigma_theta', 'sigma_Q', 'Q', 'dr_f', 'welfare', ...
    'mu_K', 'sharpe_ratio', 'excess_returns', 'nominal_interest', 'inflation', ...
    'risk_liquidity_premia', 'pdf', 'pdf_normalized', 'cdf', 'cdf_normalized', ...
    'welfare_cons',  'welfare_pi_costs'};

if s.plot_results
    % Add legends and save if needed
    figure(35)
    legend({'RP', 'LP'}, 'FontSize', s.fontsize, 'location', 'best');

    % Update leverage ylim
    figure(1);
    ylim([0, 30]);

    % Update nominal rate ylim
    figure(32);
    ylim([0, 4.5]);
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
