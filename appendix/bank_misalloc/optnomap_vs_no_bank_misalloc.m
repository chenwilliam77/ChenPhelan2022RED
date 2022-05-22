% Create plots of equilibrium dynamics when using the
% Calibration parameters
%
% Written by William Chen and Greg Phelan, Jan. 2022

close all;
addpath ../../parameters_specs
addpath(genpath('eqm_functions'))

parameters_file = 'bank_misalloc_calibration_parameters';
save_plots = 1; % Set to 1 if you want to save plots
figurespath = [pwd(), '/../../figures/bank_misalloc_optimal_fedput_nomap_vs_no_bank_misalloc/'];
load('../../save/output_data/bank_misalloc_max_expV_fedput_4args/fmincon_run_interior-point_20-Mar-2022.mat');
leg_text = {'No Misallocation', 'Misallocation'};

%% Compute equilibrium
% OptNoMaP with bank misallocation
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
[solution, summary_stats] = geteqm(s);

% Inflation targeting
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.nomrule = 2; % turn on interpolant for MP rule
[constpi_s, converged] = ...
    find_constant_inflation_mp_rule(s, 'damping', 0.7, 'tol', 5e-3);
[~, constpi_summary_stats] = geteqm(constpi_s);
s.nomrule = 1;

% OptNoMaP with no misallocation
load('../../save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat');
s = rmfield(s, 'bank_misalloc');
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
[solution2, summary_stats2] = geteqm(s);

figure(1)
plot(solution2.eta, 100 * solution2.interestvec, 'color', 'black'); hold on
plot(solution.eta, 100 * solution.interestvec, 'color', 'blue'); hold off
ylabel('\%');
ylim([0, 4.5])
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share');

figure(2)
plot(solution2.eta, 100 * solution2.inflation, 'color', 'black'); hold on
plot(solution.eta, 100 * solution.inflation, 'color', 'blue'); hold off
ylabel('annualized \%');
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share');

figure(3)
plot(solution2.eta, solution2.WelfarePi, 'color', 'black'); hold on
plot(solution.eta, solution.WelfarePi, 'color', 'blue'); hold off
ylabel('$V_\pi$');
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share');
legend(leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');

figure(4)
plot(solution2.eta, solution2.psi, 'color', 'black'); hold on
plot(solution.eta, solution.psi, 'color', 'blue'); hold off
ylabel('$\psi$');
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share');

fprintf(['The consumption equivalent percent difference b/n\n', ...
    'inflation targeting and OptNoMaP w/bank misallocation is ', ...
    num2str(round(100 * consumption_equiv(constpi_summary_stats.ergodic.ExpV, ...
    summary_stats.ergodic.ExpV, s.r), 2)), '%%\n']);


fignum = [1, 3];
save_fn = {'nominal_interest', 'WelfarePi'};

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
