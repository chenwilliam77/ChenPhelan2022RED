% Vary constant rate rules when tail risk insurance
% eliminates financial distress from equilibrium
%
% Written by William Chen, Mar. 2022

clear;
close all;
addpath('../../parameters_specs')
addpath(genpath('../../eqm_functions'))

figurespath = [pwd(), '/../../figures/equity_injection_vary_ilaw/'];
save_plots = false;

% Use equity injections to prevent financial distress from occurring
eta_TR = .06;
theta0R = 100; % make bisection faster

% Rule 1
equity_injection_calibration_parameters;
s.eta_TR = eta_TR;
s.theta0R = theta0R;
s.iput = s.ilaw;
s.lambda = 0;
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
TR1s = s;
[TR1_solution, TR1_summary_stats] = geteqm(TR1s);

% Rule 2
equity_injection_calibration_parameters;
s.eta_TR = eta_TR;
s.theta0R = theta0R;
s.iput = .02;
s.ilaw = .02;
s.lambda = 0;
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.color = 'b';
TR2s = s;
[TR2_solution, TR2_summary_stats] = geteqm(TR2s);

% Rule 3
equity_injection_calibration_parameters;
s.eta_TR = eta_TR;
s.theta0R = theta0R;
s.iput = 0;
s.ilaw = 0;
s.lambda = 0;
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.color = 'r';
TR3s = s;
[TR3_solution, TR3_summary_stats] = geteqm(TR3s);

fprintf('ExpV for i = %.1f%%:      %.4f\n', [100 * TR1s.ilaw, TR1_summary_stats.ergodic.ExpV]);
fprintf('ExpV for i = %.1f%%:     %.4f\n', [100 * TR2s.ilaw, TR2_summary_stats.ergodic.ExpV]);
fprintf('ExpV for i = %.1f%%:     %.4f\n', [100 * TR3s.ilaw, TR3_summary_stats.ergodic.ExpV]);

% Make initial plots
plot_results(TR1_solution, TR1s);
plot_results(TR2_solution, TR2s);
plot_results(TR3_solution, TR3s);

%% Edits to plots
fignum = [1, 2, 22, 23, 6, 7, 8, 9, 10, 11, 12, 16, 25, 26, 32, 33, 35, ...
    40, 41, 42, 43, 56, 58]; % figures to save/edit in loops
save_fn = {'leverage', 'mu_eta', 'sigma_eta', 'psi', 'omega', 'sigma_Ch', 'sigma_theta', 'sigma_Q', 'Q', 'dr_f', 'welfare', ...
    'mu_K', 'sharpe_ratio', 'excess_returns', 'nominal_interest', 'inflation', ...
    'risk_liquidity_premia', 'pdf', 'pdf_normalized', 'cdf', 'cdf_normalized', ...
    'welfare_cons',  'welfare_pi_costs'};

% Update nominal rate ylim
figure(32);
ylim([0, 5]);

% Add legends and save if needed
for idx = 1:length(fignum)
    figure(fignum(idx));
    legend({'No Tail Risk', ...
        ['$i = ' num2str(round(TR2s.ilaw * 100, 0)) '$\%'], ...
        ['$i = ' num2str(round(TR3s.ilaw * 100, 0)) '$\%']}, 'FontSize', s.fontsize, 'location', 'best');
end

if save_plots
    if ~isfolder(figurespath)
        mkdir(figurespath);
    end

    for i = 1:length(fignum)
        f = figure(fignum(i));
        saveas(f, [figurespath, '/', save_fn{i}, '.', s.graphics_format]);
    end
end
