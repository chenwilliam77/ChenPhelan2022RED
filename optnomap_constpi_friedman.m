% Create plots of equilibrium dynamics under three policies:
% (1) Fed Put
% (2) constant inflation
% (3) Friedman Rule adjusted for inflation costs in first-best eqm
%
% Written by William Chen and Gregory Phelan, Feb. 2022

clear;
close all;
addpath parameters_specs
addpath eqm_functions

parameters_file = 'calibration_parameters';
figurespath = [pwd(), '/figures/optnomap_constpi_friedman/'];
save_plots = 1; % Set to 1 if you want to save plots
find_constpi_rule = 1; % Set to 1 to calculate the constant inflation MP rule, otherwise load from memory
constpi_path = [pwd(), '/save/output_data/fedput_calibration_constpi/']; % save path for constant inflation MP rule
load(['save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat']);

leg_text = {'OptNoMaP', 'Inflation Targeting', 'Friedman, $\lambda > 0$'};
lvg_leg_text = {'Inflation Targeting', '$\hat{i}^{Put}$', '$\hat{i}^{LAW}$'};

%% Compute equilibrium
% Fed put
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
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

% Friedman taking into account inflation costs
eval(parameters_file);
constrate = compute_first_best_rate(s);
s.ilaw = constrate;
s.ireg = constrate;
s.iput = constrate;
s.plot_results = 0;
s.use_title = 0;
s.color = 'r';
friedman_s = s;
[friedman_solution, friedman_summary_stats] = geteqm(friedman_s);

% Make initial plots
plot_results(fp_solution, fp_s);
plot_results(constpi_solution, constpi_s);
plot_results(friedman_solution, friedman_s);

%% Edits to plots
fignum = [1, 2, 22, 23, 6, 7, 8, 9, 10, 11, 12, 16, 25, 26, 32, 33, 35, ...
    40, 41, 42, 43, 56, 58]; % figures to save/edit in loops
save_fn = {'leverage', 'mu_eta', 'sigma_eta', 'psi', 'omega', 'sigma_Ch', 'sigma_theta', 'sigma_Q', 'Q', 'dr_f', 'welfare', ...
    'mu_K', 'sharpe_ratio', 'excess_returns', 'nominal_interest', 'inflation', ...
    'risk_liquidity_premia', 'pdf', 'pdf_normalized', 'cdf', 'cdf_normalized', ...
    'welfare_cons',  'welfare_pi_costs'};

% Add legends and save if needed
for idx = 1:length(fignum)
    figure(fignum(idx));
    if fignum(idx) == 35
        legend({'HP + SP', 'LP'}, 'FontSize', s.legend_fontsize, 'location', 'best');
    else
        legend(leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');
    end
end

% Update leverage ylim
figure(1);
ylim([0, 30]);

% Update nominal rate ylim
figure(32);
ylim([0, 4.5]);

% Plot leverage as percentage differences from OptNoMaP
% First compute ilaw and iput constant cases
eval(parameters_file);
s.ilaw = best_rule(1);
s.iput = best_rule(1);
s.plot_results = 0;
s.use_title = 0;
iput_s = s;
[iput_solution, iput_summary_stats] = geteqm(iput_s);

s.ilaw = best_rule(2);
s.iput = best_rule(2);
s.plot_results = 0;
s.use_title = 0;
ilaw_s = s;
[ilaw_solution, ilaw_summary_stats] = geteqm(ilaw_s);

figure(900);
interp_fp_lvg_friedman = interp1(fp_solution.eta, fp_solution.Leverage, friedman_solution.eta, 'linear');
interp_fp_lvg_iput = interp1(fp_solution.eta, fp_solution.Leverage, iput_solution.eta, 'linear');
interp_fp_lvg_ilaw = interp1(fp_solution.eta, fp_solution.Leverage, ilaw_solution.eta, 'linear');
interp_fp_lvg_constpi = interp1(fp_solution.eta, fp_solution.Leverage, constpi_solution.eta, 'linear');
plot(constpi_solution.eta, 100 * (constpi_solution.Leverage ./ interp_fp_lvg_constpi - 1), 'color', constpi_s.color, 'LineWidth', fp_s.width); hold on
plot(iput_solution.eta, 100 * (iput_solution.Leverage ./ interp_fp_lvg_iput - 1), 'color', [0.8500 0.3250 0.0980], 'LineWidth', fp_s.width); hold on
plot(ilaw_solution.eta, 100 * (ilaw_solution.Leverage ./ interp_fp_lvg_ilaw - 1), 'color', [0.4940 0.1840 0.5560], 'LineWidth', fp_s.width); hold on
max_eta_lvg_plot = max([fp_solution.eta(end), constpi_solution.eta(end), ...
    iput_solution.eta(end), ilaw_solution.eta(end)]);
plot([0, max_eta_lvg_plot], [0, 0], 'color', 'black', 'LineStyle', '--');
xlim([0, max_eta_lvg_plot]);
legend(lvg_leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');
ylabel('\% Difference');
xlabel('$\eta$')

% Save plots
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

% Plot consumption equivalent differences
welfare_summary_stats = compare_welfare({fp_solution, constpi_solution, friedman_solution}, ...
                                        {fp_summary_stats, constpi_summary_stats, friedman_summary_stats}, ...
                                        1000, leg_text(2:end), ['b', 'r'], ...
                                        fp_s, save_plots, figurespath);
