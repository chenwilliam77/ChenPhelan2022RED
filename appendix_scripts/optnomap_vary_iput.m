% Create plots of equilibrium dynamics under two policies:
% (1) Optimal Fed Put (No MaP)
% (2) iput higher than Optimal Fed Put
%
% Written by William Chen and Gregory Phelan, Mar. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

figurespath = [pwd(), '/../figures/optnomap_vary_iput/'];
save_plots = 1; % Set to 1 if you want to save plots

%% Compute equilibrium
% Optimal Fed Put (No MaP)
load('../save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat');
calibration_parameters;
s.plot_results = 0;
s.use_title = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
fp_s = s; % copy the struct
[fp_solution, fp_summary_stats] = geteqm(fp_s);

% Increase iput
s.iput = best_rule(1) + .03;
s.color = 'b';
mppol2_s = s;
[mppol2_solution, mppol2_summary_stats] = geteqm(mppol2_s);

% Make initial plots
plot_results(fp_solution, fp_s);
plot_results(mppol2_solution, mppol2_s);

%% Edits to plots
fignum = [26, 35]; % figures to save/edit in loops
save_fn = {'excess_returns', 'risk_liquidity_premia'};

% Figure out what the second policy's nominal rate should be printed as
sec_pol_txt = 'Higher $i^{Put}$';

% Add legends and save if needed
for idx = 1:length(fignum)
    figure(fignum(idx));
    if fignum(idx) == 35
        legend({'HP + SP', 'LP'}, 'FontSize', s.legend_fontsize, 'location', 'best');
    else
        legend({'OptNoMAP', sec_pol_txt}, 'FontSize', s.legend_fontsize, 'location', 'best');
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
