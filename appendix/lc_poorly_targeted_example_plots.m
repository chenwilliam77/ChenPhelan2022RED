% Plot leverage under a poorly targeted
% macroprudential policy binding during
% the Crisis region
%
% Written by William Chen and Gregory Phelan, Feb. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

parameters_file = 'lc_linear_calibration_parameters';
do_plot_results = 1;
figurespath = [pwd(), '/../figures/lc_linear_crisis_optimal_fedput_nomap/'];
save_plots = 1;
if save_plots && ~isfolder(figurespath)
    mkdir(figurespath);
end
load('../save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat');
leg_text = {'No MaP', 'Looser MaP', 'Tighter MaP'};

% Policy 1
eval(parameters_file);
s.leverage_constraint.type = 'linear';
s.leverage_constraint.etagrid = [0; .02; 1];
s.leverage_constraint.Ls = [Inf; Inf; Inf];
s.plot_results = 0;
s.use_title = 0;
s.autodiff_jacobian = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol1s = s;
[pol1_solution, polInf_summary_stats] = geteqm(pol1s);

% Policy 2
eval(parameters_file);
s.leverage_constraint.type = 'linear';
s.leverage_constraint.etagrid = [0; .02; 1];
s.leverage_constraint.Ls = [35; 25; Inf];
s.plot_results = 0;
s.use_title = 0;
s.color = 'b';
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol2s = s;
[pol2_solution, pol2_summary_stats] = geteqm(pol2s);

% Policy 3
eval(parameters_file);
s.leverage_constraint.type = 'linear';
s.leverage_constraint.etagrid = [0; .02; 1];
s.leverage_constraint.Ls = [30; 25; Inf];
s.plot_results = 0;
s.use_title = 0;
s.color = 'r';
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol3s = s;
[pol3_solution, pol3_summary_stats] = geteqm(pol3s);

fprintf('ExpV for L = %.1f:      %.2f\n', [pol1s.leverage_constraint.Ls(2), polInf_summary_stats.ergodic.ExpV]);
fprintf('ExpV for L = %.1f:      %.2f\n', [pol2s.leverage_constraint.Ls(2), pol2_summary_stats.ergodic.ExpV]);
fprintf('ExpV for L = %.1f:      %.2f\n', [pol3s.leverage_constraint.Ls(2), pol3_summary_stats.ergodic.ExpV]);
fprintf('ExpVc for L = %.1f:     %.2f\n', [pol1s.leverage_constraint.Ls(2), polInf_summary_stats.ergodic.ExpVc]);
fprintf('ExpVc for L = %.1f:     %.2f\n', [pol2s.leverage_constraint.Ls(2), pol2_summary_stats.ergodic.ExpVc]);
fprintf('ExpVc for L = %.1f:     %.2f\n', [pol3s.leverage_constraint.Ls(2), pol3_summary_stats.ergodic.ExpVc]);

if do_plot_results
% Make initial plots
plot_results(pol1_solution, pol1s);
plot_results(pol2_solution, pol2s);
plot_results(pol3_solution, pol3s);

%% Edits to plots
fignum = [1]; % figures to edit in loops
save_fn = {'leverage'};

figure(1);
ylim([0, 60]);

% Add legends and save if needed
for idx = 1:length(fignum)
    figure(fignum(idx));
    legend(leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');
end

if save_plots
    for i = 1:length(fignum)
        f = figure(fignum(i));
        saveas(f, [figurespath, '/', save_fn{i}, '.', s.graphics_format]);
    end
end
end
