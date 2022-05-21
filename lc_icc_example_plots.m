% Plot the stationary density and leverage percentage difference
% under different ICC rules while holding fixed MP
% at OptNoMaP.
%
% Written by William Chen and Gregory Phelan, Mar. 2022

clear;
close all;
addpath parameters_specs
addpath(genpath('eqm_functions'))

parameters_file = 'lc_icc_calibration_parameters';
do_plot_results = 1;
figurespath = [pwd(), '/figures/lc_icc_example_plots/'];
save_plots = 1;
load(['save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat']);
leg_text = {'$\alpha = \infty$', '$\alpha = 10$', '$\alpha = 8$'};

% Policy 1
eval(parameters_file);
s.leverage_constraint.alpha = Inf;
s.plot_results = 0;
s.use_title = 0;
s.autodiff_jacobian = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol1s = s;
[pol1_solution, pol1_summary_stats] = geteqm(pol1s);

% Policy 2
eval(parameters_file);
s.leverage_constraint.type = 'icc';
s.leverage_constraint.alpha = 10;
s.plot_results = 0;
s.use_title = 0;
s.color = 'b';
s.autodiff_jacobian = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol2s = s;
[pol2_solution, pol2_summary_stats] = geteqm(pol2s);

% Policy 3
eval(parameters_file);
s.leverage_constraint.type = 'icc';
s.leverage_constraint.alpha = 8;
s.plot_results = 0;
s.use_title = 0;
s.color = 'r';
s.autodiff_jacobian = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol3s = s;
[pol3_solution, pol3_summary_stats] = geteqm(pol3s);

fprintf('ExpV for alpha = %.1f:      %.2f\n', [pol1s.leverage_constraint.alpha, pol1_summary_stats.ergodic.ExpV]);
fprintf('ExpV for alpha = %.1f:      %.2f\n', [pol2s.leverage_constraint.alpha, pol2_summary_stats.ergodic.ExpV]);
fprintf('ExpV for alpha = %.1f:      %.2f\n', [pol3s.leverage_constraint.alpha, pol3_summary_stats.ergodic.ExpV]);
fprintf('ExpVc for alpha = %.1f:     %.2f\n', [pol1s.leverage_constraint.alpha, pol1_summary_stats.ergodic.ExpVc]);
fprintf('ExpVc for alpha = %.1f:     %.2f\n', [pol2s.leverage_constraint.alpha, pol2_summary_stats.ergodic.ExpVc]);
fprintf('ExpVc for alpha = %.1f:     %.2f\n', [pol3s.leverage_constraint.alpha, pol3_summary_stats.ergodic.ExpVc]);

if do_plot_results
% Make initial plots
plot_results(pol1_solution, pol1s);
plot_results(pol2_solution, pol2s);
plot_results(pol3_solution, pol3s);

%% Edits to plots
fignum = [40, 41, 42, 43]; % figures to edit in loops
save_fn = {'pdf', 'pdf_normalized', 'cdf', 'cdf_normalized'};

% Add legends and save if needed
for idx = 1:length(fignum)
    figure(fignum(idx));
    legend(leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');
end

% Plot leverage as percentage differences from baseline policy
figure(900);
interp_iccInf_lvg_icc1p5 = interp1(pol1_solution.eta, pol1_solution.Leverage, pol2_solution.eta, 'linear');
interp_iccInf_lvg_icc1p2 = interp1(pol1_solution.eta, pol1_solution.Leverage, pol3_solution.eta, 'linear');
plot(pol2_solution.eta(10:end), 100 * (pol2_solution.Leverage(10:end) ./ interp_iccInf_lvg_icc1p5(10:end) - 1), 'color', pol2s.color, 'LineWidth', s.width); hold on
plot(pol3_solution.eta(10:end), 100 * (pol3_solution.Leverage(10:end) ./ interp_iccInf_lvg_icc1p2(10:end) - 1), 'color', pol3s.color, 'LineWidth', s.width); hold on
max_eta = max([pol1_solution.eta(end), pol2_solution.eta(end), pol3_solution.eta(end)]);
plot([0, max_eta], [0, 0], 'color', 'black', 'LineStyle', '--');
legend(leg_text(2:end), 'FontSize', s.fontsize');
xlim([0, max_eta]);
ylabel('\% Difference');
xlabel('$\eta$')

if save_plots
    if ~isfolder(figurespath)
        mkdir(figurespath);
    end
    for i = 1:length(fignum)
        f = figure(fignum(i));
        saveas(f, [figurespath, '/', save_fn{i}, '.', s.graphics_format]);
    end
    f = figure(900);
    saveas(f, [figurespath, '/leverage_pctdev.', s.graphics_format]);
end
end
