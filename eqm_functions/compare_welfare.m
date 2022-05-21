function welfare_summary_stats = compare_welfare(sol_cell, sum_stat_cell, first_fignum, legend_text, colors, s, save_plots, figurespath);
% welfare_summary_stats = compare_welfare(sol_cell, sum_stat_cell, first_fignum, legend_text, colors, s, save_plots, figurespath);
%
% Compares welfare across different equilibria through
% plots and summary statistics
%
% sol_cell:      cell of solution structs; it is assumed
%                that the first struct is the reference against
%                which comparisons are made
%
% sum_stat_cell: cell of summary statistics structs; it is assumed
%                that the first struct is the reference against
%                which comparisons are made
%
% first_fignum:  starting figure number for plots
%
% legend_text:   text for the legend
%
% colors:        colors for the lines
%
% s:             settings struct
%
% save_plots:    indicator for whether to save plots (optional)
%
% figurespath:   path to the figures folder to save plots (optional)
%
% Written by William Chen and Gregory Phelan, Jan. 2022

% Set up
N = length(sol_cell) - 1; % number of things to plots, not number of equilibria
etavec = sol_cell{1}.eta;
V_ref_eta_cell = cell(N, 1);  % welfare of reference equilibrium evaluated at eta
Vc_ref_eta_cell = cell(N, 1); % for the other equilibria.
Vd_ref_eta_cell = cell(N, 1); % Each cell's element holds a vector
Vpi_ref_eta_cell = cell(N, 1);
V_comp_cell = cell(N, 1);  % welfare comparisons of other equilibria against
Vc_comp_cell = cell(N, 1); % the reference equilibrium for various objects
Vc_Vpi_comp_cell = cell(N, 1);
other_eqm_max_eta_inds = zeros(N, 1); % index of last eta in other equilibria that also belongs in etavec

% Project welfare from reference equilibrium onto state vector
% for the other equilibria using interpolation
for i = 1:N
    if etavec(end) >= sol_cell{i + 1}.eta(end) % in this case, we can directly interpolate onto other eqm's eta
        other_eqm_max_eta_inds(i) = length(sol_cell{i + 1}.eta);
        V_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.Welfare, sol_cell{i + 1}.eta, 'linear');
        Vc_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.WelfareConsume, sol_cell{i + 1}.eta, 'linear');
        Vd_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.WelfareDeposits, sol_cell{i + 1}.eta, 'linear');
        Vpi_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.WelfarePi, sol_cell{i + 1}.eta, 'linear');
    else % other eqm has longer eta => need to truncate evaluation
        other_eqm_max_eta_inds(i) = find(sol_cell{i + 1}.eta <= etavec(end), 1, 'last');
        V_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.Welfare, ...
            sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), 'linear');
        Vc_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.WelfareConsume, ...
            sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), 'linear');
        Vd_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.WelfareDeposits, ...
            sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), 'linear');
        Vpi_ref_eta_cell{i} = interp1(etavec, sol_cell{1}.WelfarePi, ...
            sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), 'linear');
    end
end

% Calculate consumption equivalents in percentage terms
for i = 1:N
    other_max_eta_ind = other_eqm_max_eta_inds(i);
    for j = 1:other_max_eta_ind
        V_comp_cell{i}(j) = 100 * consumption_equiv(sol_cell{i + 1}.Welfare(j), V_ref_eta_cell{i}(j), s.r);
        Vc_comp_cell{i}(j) = 100 * consumption_equiv(sol_cell{i + 1}.WelfareConsume(j), Vc_ref_eta_cell{i}(j), s.r);
        Vc_Vpi_comp_cell{i}(j) = 100 * consumption_equiv(sol_cell{i + 1}.Welfare(j) - sol_cell{i + 1}.WelfareDeposits(j), ...
            V_ref_eta_cell{i}(j) - Vd_ref_eta_cell{i}(j), s.r);
    end
end

% Plot consumption equivalent differences

% V
figure(first_fignum); hold on
for i = 1:N
    plot(sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), V_comp_cell{i}, 'color', colors(i), 'linestyle', s.line, 'linewidth', s.width);
end
legend(legend_text, 'FontSize', s.fontsize, 'location', 'best'); hold off
eta_ends = zeros(N, 1);
min_V_vec = zeros(N, 1);
max_V_vec = zeros(N, 1);
for i = 1:N
    eta_ends(i) = sol_cell{i}.eta(end);
    min_V_vec(i) = min(V_comp_cell{i});
    max_V_vec(i) = max(V_comp_cell{i});
end
xlim([s.start; 1.05 * max(eta_ends)]);
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylabel('\% Difference', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylim(1.05 * [min(min(min_V_vec)),  max([max_V_vec; 0])]);

% Vc
figure(first_fignum + 1); hold on
for i = 1:N
    plot(sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), Vc_comp_cell{i}, 'color', colors(i), 'linestyle', s.line, 'linewidth', s.width);
end
legend(legend_text, 'FontSize', s.fontsize, 'location', 'best'); hold off
min_Vc_vec = zeros(N, 1);
max_Vc_vec = zeros(N, 1);
for i = 1:N
    min_Vc_vec(i) = min(Vc_comp_cell{i});
    max_Vc_vec(i) = max(Vc_comp_cell{i});
end
xlim([s.start; 1.05 * max(eta_ends)]);
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylabel('\% Difference', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylim(1.05 * [min(min(min_Vc_vec)),  max([max_Vc_vec; 0])]);

% Vc - Vpi
figure(first_fignum + 2); hold on
for i = 1:N
    plot(sol_cell{i + 1}.eta(1:other_eqm_max_eta_inds(i)), Vc_Vpi_comp_cell{i}, 'color', colors(i), 'linestyle', s.line, 'linewidth', s.width);
end
legend(legend_text, 'FontSize', s.fontsize, 'location', 'best'); hold off
min_Vc_Vpi_vec = zeros(N, 1);
max_Vc_Vpi_vec = zeros(N, 1);
for i = 1:N
    min_Vc_Vpi_vec(i) = min(Vc_Vpi_comp_cell{i});
    max_Vc_Vpi_vec(i) = max(Vc_Vpi_comp_cell{i});
end
xlim([s.start; 1.05 * max(eta_ends)]);
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylabel('\% Difference', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylim(1.05 * [min(min_Vc_Vpi_vec); max([max_Vc_Vpi_vec; 0])]);

% Calculate summary statistics
expV_vec = zeros(N, 1);
expVc_vec = zeros(N, 1);
expVc_Vpi_vec = zeros(N, 1);
for i = 1:N
    expV_vec(i) = 100 * consumption_equiv(sum_stat_cell{i + 1}.ergodic.ExpV, sum_stat_cell{1}.ergodic.ExpV, s.r);
    expVc_vec(i) = 100 * consumption_equiv(sum_stat_cell{i + 1}.ergodic.ExpVc, sum_stat_cell{1}.ergodic.ExpVc, s.r);
    expVc_Vpi_vec(i) = 100 * consumption_equiv(sum_stat_cell{i + 1}.ergodic.ExpV - sum_stat_cell{i + 1}.ergodic.ExpVd, ...
        sum_stat_cell{1}.ergodic.ExpV - sum_stat_cell{1}.ergodic.ExpVd, s.r);
end

% Save plots if desired
if nargin > 6
    if save_plots
	    f = figure(first_fignum);
	    saveas(f, [figurespath, '/consumption_equiv_V.', s.graphics_format]);
	    f = figure(first_fignum + 1);
	    saveas(f, [figurespath, '/consumption_equiv_Vc.', s.graphics_format]);
	    f = figure(first_fignum + 2);
	    saveas(f, [figurespath, '/consumption_equiv_Vc_Vpi.', s.graphics_format]);
    end
end

% Return output struct
welfare_summary_stats.expV_vec = expV_vec;
welfare_summary_stats.expVc_vec = expVc_vec;
welfare_summary_stats.expVc_Vpi_vec = expVc_Vpi_vec;

end
