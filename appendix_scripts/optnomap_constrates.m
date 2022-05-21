% Plot OptNoMaP vs. constant rate rules (\tilde{i}^{Put}, \tilde{i}^{LAW})
%
% Written by William Chen and Gregory Phelan, Mar. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

parameters_file = 'calibration_parameters';
load('../save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat');

% Rule 1 - Optimal Fed Put
n_comps = 5;
rule1_stats.welfare = zeros(n_comps, 1);
rule1_stats.stab_frac = zeros(n_comps, 1);
rule1_stats.stab_frac_pct50 = zeros(n_comps, 1);
rule1_sols = cell(n_comps, 1);
rule1_sumstats = cell(n_comps, 1);

eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;

% Best rule
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);

L1s = s;
[rule1_sols{1}, rule1_sumstats{1}] = geteqm(L1s);

% Make initial plots
plot_results(rule1_sols{1}, L1s);


% OptNoMaP
rule1_stats.welfare(1) = 0;
rule1_stats.stab_frac(1) = 0;

% i = iput constant
L1s.ilaw = best_rule(1);
[rule1_sols{2}, rule1_sumstats{2}] = geteqm(L1s);

% Make initial plots
L1s.color='b';
plot_results(rule1_sols{2}, L1s);

% i = ilaw constant
L1s.ilaw = best_rule(2);
L1s.iput = best_rule(2);
[rule1_sols{3}, rule1_sumstats{3}] = geteqm(L1s);

% Make initial plots
L1s.color='r';
plot_results(rule1_sols{3}, L1s);

%%
leg_text = {'OptNoMaP', '$\tilde{i}^{Put}$', '$\tilde{i}^{LAW}$'};

figurespath = [pwd(), '/../figures/optnomap_constrates/'];
if ~isfolder(figurespath)
    mkdir(figurespath);
end
fignum = [9,2,22,35];
save_fn = {'sigma_Q_constrates','mu_eta_constrates','sigma_eta_constrates','risk_liquidity_premia_constrates'};


for i=1:length(fignum)-1
    f = figure(fignum(i));
    if fignum(i) == 9
        legend(leg_text, 'FontSize', s.legend_fontsize, 'location', 'best');
    end
    saveas(f, [figurespath, '/', save_fn{i}, '.', s.graphics_format]);
end

f = figure(fignum(length(fignum)));
saveas(f, [figurespath, '/', save_fn{length(fignum)}, '.', s.graphics_format]);
