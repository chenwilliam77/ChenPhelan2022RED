% Vary alpha given OptMP, inflation targeting, and Friedman (lambda > 0) rules
%
% Written by William Chen and Gregory Phelan, Mar. 2022

clear;
close all;
addpath parameters_specs
addpath(genpath('eqm_functions'))

tic

parameters_file = 'lc_icc_calibration_parameters';
do_plot_results = 1;
figurespath = [pwd(), '/figures/lc_icc_optmp_vary_alpha/'];
save_plots = 1;
leg_text = {'OptMP', 'Inflation Targeting', 'Friedman, $\lambda > 0$'};
colors = {'black', 'blue', 'red'};

alphavec = [14, 12, 10, 8, 6, 4];
day_cells = {''; ''; '08'; '08'; '08'; '12'};
CE = alphavec;
distress = CE;
crisis = CE;

% Policy 1
eval(parameters_file);
load('save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat');
s.leverage_constraint.alpha = Inf;
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol1s = s;
[pol_solution, pol_summary_stats] = geteqm(pol1s);

friedman_rate_lambda_pos = compute_first_best_rate(pol1s);
pol1s.ilaw = friedman_rate_lambda_pos;
pol1s.iput = friedman_rate_lambda_pos;

% Inflation Targeting, save for later
[constpi_s, converged] = ...
    find_constant_inflation_mp_rule(pol1s, 'damping', 0.7, 'tol', 5e-3);
if ~converged
    disp('Inflation targeting rule not found! Figure out what is wrong');
    keyboard;
end

% Now loop over alpha
eval(parameters_file);
s.plot_results = 0;
s.use_title = 0;
s.verbose = 0;
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
s.leverage_constraint.type = 'icc';

% first we run with MP=OptMP
pol_alpha_s = cell(length(alphavec), 1);
pol_alpha_sol = cell(length(alphavec), 1);
pol_alpha_sumstats = cell(length(alphavec), 1);
for i=1:length(alphavec)
if i == 1 || i == 2 % first two alphas are "dummies" that shouldn't bind
    load(['save/output_data/max_expV_fedput_4args/fmincon_run_interior-point_07-Mar-2022.mat']);
else
	load(['save/output_data/max_expV_fedput_4args_icc' num2str(alphavec(i)) ...
	    '/fmincon_run_interior-point_' day_cells{i} '-Mar-2022.mat']);
end

% Policy w alpha
s.iput = best_rule(1);
s.ilaw = best_rule(2);
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol_alpha_s{i} = s;
pol_alpha_s{i}.leverage_constraint.alpha = alphavec(i);

[pol_alpha_sol{i}, pol1_alpha_sumstats{i}] = geteqm(pol_alpha_s{i});

CE(i)=100*consumption_equiv( pol1_alpha_sumstats{i}.ergodic.ExpV, pol_summary_stats.ergodic.ExpV, s.r);
distress(i)=100 - pol1_alpha_sumstats{i}.ergodic.stab_frac;
crisis(i)=100 - pol1_alpha_sumstats{i}.ergodic.stab_frac_pct50;
end

figure(1)
plot(alphavec,CE,'-x', 'MarkerSize', 10, 'MarkerFaceColor', colors{1}, ...
    'Color', colors{1}); hold on
ylabel('Welfare Gains')
xlabel('$\alpha$, tightness of ICC constraint')

figure(2)
plot(alphavec,distress, '-x', 'MarkerSize', 10, 'MarkerFaceColor', colors{1}, ...
    'Color', colors{1}); hold on
ylabel('Probability of Distress')
xlabel('$\alpha$, tightness of ICC constraint')

figure(3)
plot(alphavec,crisis, '-x', 'MarkerSize', 10, 'MarkerFaceColor', colors{1}, ...
    'Color', colors{1}); hold on
ylabel('Probability of Crisis')
xlabel('$\alpha$, tightness of ICC constraint')

toc

%
% second we run with 'inflation targeting'

CE_IT=NaN*CE;
distress_IT=NaN*CE;
crisis_IT=NaN*CE;

for i=1:length(alphavec)
% Policy w alpha
if i==1
s=pol1s(i);
else
    s=constpi_s;
end
s.ilaw = friedman_rate_lambda_pos;
s.iput = friedman_rate_lambda_pos;
s.leverage_constraint.type = 'icc';
s.leverage_constraint.alpha = alphavec(i);
if alphavec(i) < 11 && alphavec(i) ~= 6
    tol=1e-2;
    thinning = 5;
    damping = 0.7;
elseif alphavec(i) == 6
    tol = 1e-2;
    thinning = 5;
    damping = 0.8;
else
    tol=5e-3;
    thinning = 10;
    damping = 0.7;
end

% Inflation Targeting, save for later
[constpi_s, converged] = ...
    find_constant_inflation_mp_rule(s, 'damping', damping, 'tol', tol, 'thinning', thinning);
disp(['Inflation targeting for alpha = ', ...
    num2str(alphavec(i)), ' has convergence code ' num2str(converged)]);

if ~converged
    keyboard;
end
constpi_s.leverage_constraint.type = 'icc';
constpi_s.leverage_constraint.alpha = alphavec(i);
[pol_alpha_sol{i}, pol1_alpha_sumstats{i}] = geteqm(constpi_s);

CE_IT(i)=100*consumption_equiv( pol1_alpha_sumstats{i}.ergodic.ExpV, pol_summary_stats.ergodic.ExpV, s.r);
distress_IT(i)=100 - pol1_alpha_sumstats{i}.ergodic.stab_frac;
crisis_IT(i)=100 - pol1_alpha_sumstats{i}.ergodic.stab_frac_pct50;
end

figure(1)
plot(alphavec,CE_IT, '-h', 'MarkerSize', 7, 'MarkerFaceColor', colors{2}, ...
    'Color', colors{2}); hold on
ylabel('Welfare Gains')
xlabel('$\alpha$, tightness of ICC constraint')

figure(2)
plot(alphavec,distress_IT, '-h', 'MarkerSize', 7, 'MarkerFaceColor', colors{2}, ...
    'Color', colors{2}); hold on
ylabel('Probability of Distress')
xlabel('$\alpha$, tightness of ICC constraint')

figure(3)
plot(alphavec,crisis_IT, '-h', 'MarkerSize', 7, 'MarkerFaceColor', colors{2}, ...
    'Color', colors{2}); hold on
ylabel('Probability of Crisis')
xlabel('$\alpha$, tightness of ICC constraint')

toc

%
%third we run with MP=Friedman
eval(parameters_file);
friedman_rate_lambda_pos = compute_first_best_rate(s);

s.leverage_constraint.type = 'icc';
s.plot_results = 0;
s.use_title = 0;
s.autodiff_jacobian = 0;
s.verbose = 0;
s.iput = friedman_rate_lambda_pos;
s.ilaw = friedman_rate_lambda_pos;
s.etaPUT = best_rule(3);
s.etaLAW = best_rule(4) + best_rule(3);
pol1s = s;

CE_F=NaN*CE;
distress_F=NaN*CE;
crisis_F=NaN*CE;

for i=1:length(alphavec)

% Policy w alpha
pol_alpha_s{i} = s;
pol_alpha_s{i}.leverage_constraint.alpha = alphavec(i);

[pol_alpha_sol{i}, pol1_alpha_sumstats{i}] = geteqm(pol_alpha_s{i});

CE_F(i)=100*consumption_equiv( pol1_alpha_sumstats{i}.ergodic.ExpV, pol_summary_stats.ergodic.ExpV, s.r);
distress_F(i)=100 - pol1_alpha_sumstats{i}.ergodic.stab_frac;
crisis_F(i)=100 - pol1_alpha_sumstats{i}.ergodic.stab_frac_pct50;

end

figure(1)
plot(alphavec,CE_F, '-|', 'MarkerSize', 10, 'MarkerFaceColor', colors{3}, ...
    'Color', colors{3}); hold on
ylabel('Welfare Gains')
xlabel('$\alpha$, tightness of ICC constraint')
legend(leg_text, 'location', 'best');

figure(2)
plot(alphavec,distress_F, '-|', 'MarkerSize', 10, 'MarkerFaceColor', colors{3}, ...
    'Color', colors{3}); hold on
ylabel('Probability of Distress')
xlabel('$\alpha$, tightness of ICC constraint')

figure(3)
plot(alphavec,crisis_F, '-|', 'MarkerSize', 10, 'MarkerFaceColor', colors{3}, ...
    'Color', colors{3}); hold on
ylabel('Probability of Crisis')
xlabel('$\alpha$, tightness of ICC constraint')

toc

%%
if save_plots
    if ~isfolder(figurespath)
        mkdir(figurespath);
    end
    f = figure(1);
    saveas(f, [figurespath, '/Welfare_vary_alpha.', s.graphics_format]);

    f = figure(2);
    saveas(f, [figurespath, '/Distress_vary_alpha.', s.graphics_format]);

    f = figure(3);
    saveas(f, [figurespath, '/Crisis_vary_alpha.', s.graphics_format]);
end
