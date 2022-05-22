% Create plots of equilibrium dynamics when using the
% Calibration parameters
%
% Written by William Chen and Greg Phelan, Jan. 2022

clear;
close all;
addpath ../../parameters_specs
addpath(genpath('eqm_functions'))

parameters_file = 'bank_misalloc_calibration_parameters';
save_plots = 1; % Set to 1 if you want to save plots
figurespath = [pwd(), '/../../figures/bank_misalloc_fedput_calibration/'];

%% Compute equilibrium
% Fed put
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.crisis_output_losses = 0; % report crisis output losses in the summary statistics
rng(1793);
[solution, summary_stats] = geteqm(s);
s.crisis_output_losses = 0;

[max_gdp_val, max_gdp] = max(solution.gdp);
min_gdp_val = min(solution.gdp);
min_optim_psi = find(solution.gdp >= solution.gdp(end), 1);

figure(1)
yyaxis left
plot(solution.eta, solution.gdp, 'color', 'black');
patch('XData', [solution.eta(min_optim_psi), solution.eta(max_gdp), ...
    solution.eta(max_gdp), solution.eta(min_optim_psi)], ...
    'YData', [min_gdp_val, min_gdp_val, max_gdp_val, max_gdp_val], ...
    'FaceColor', 'black', 'FaceAlpha', .3)
ylabel('Output', 'color', 'black'); hold on
hold off
ylim([min_gdp_val, max_gdp_val])
yyaxis right
plot(solution.eta, solution.psi, 'color', 'black', 'LineStyle', '--');
ylh = ylabel('$\psi$', 'color', 'black');
set(ylh, 'rotation', 270);
ylim([0, 1])
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

fignum = [1];
save_fn = {'gdp'};

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
