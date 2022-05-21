function stab_summary_stats = compare_stability(sol_cell, sum_stat_cell, first_fignum, legend_text, colors, s, save_plots, figurespath);
% stab_summary_stats = compare_stability(sol_cell, sum_stat_cell, first_fignum, legend_text, colors, s, save_plots, figurespath);
%
% Compares stability across different equilibria through
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
% Written by William Chen and Gregory Phelan, Feb. 2022

% Set up
N = length(sol_cell) - 1; % number of things to plots, not number of equilibria
cdf_interps = cell(N, 1);
pdf_interps = cell(N, 1);
cdf_diff = cell(N, 1);
pdf_diff = cell(N, 1);

% Create interpolation functions
for i = 1:N
    sol = sol_cell{i + 1};
    eta_normalized = sol.eta_density ./ sol.eta_density(end);
    cdf_interps{i} = griddedInterpolant(eta_normalized, sol.cdf);
    pdf_interps{i} = griddedInterpolant(eta_normalized, sol.eta_density(end) * sol.density);
end

% Compute relative differences of PDF/CDF
sol1 = sol_cell{1};
xvec = sol1.eta_density ./ sol1.eta_density(end);
pdf_norm1 = sol1.density * sol1.eta_density(end);
for i = 1:N
    cdf_diff{i} = cdf_interps{i}(xvec) - sol1.cdf;
    pdf_diff{i} = pdf_interps{i}(xvec) - pdf_norm1;
end

% Plot relative differences

% PDF
figure(first_fignum); hold on
for i = 1:N
    plot(xvec, pdf_diff{i}, 'color', colors(i), 'linestyle', s.line, 'linewidth', s.width);
end
legend(legend_text, 'FontSize', s.fontsize, 'location', 'best'); hold off
min_pdf_vec = zeros(N, 1);
max_pdf_vec = zeros(N, 1);
for i = 1:N
    min_pdf_vec(i) = min(pdf_diff{i});
    max_pdf_vec(i) = max(pdf_diff{i});
end
xlim([0, 1]);
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylabel('Difference in $f(\eta / \eta^*)$', 'FontSize', s.fontsize, 'FontWeight', 'b');
ylim([min([0, 1.05 * min(min_pdf_vec)]), max([0, 1.05 * max(max_pdf_vec)])]);

% CDF
figure(first_fignum + 1); hold on
for i = 1:N
    plot(xvec, cdf_diff{i}, 'color', colors(i), 'linestyle', s.line, 'linewidth', s.width);
end
legend(legend_text, 'FontSize', s.fontsize, 'location', 'best'); hold off
min_cdf_vec = zeros(N, 1);
max_cdf_vec = zeros(N, 1);
for i = 1:N
    min_cdf_vec(i) = min(cdf_diff{i});
    max_cdf_vec(i) = max(cdf_diff{i});
end
xlim([0, 1]);
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylabel('Difference in $F(\eta / \eta^*)$', 'FontSize', s.fontsize, 'FontWeight', 'b');
ylim([min([0, 1.05 * min(min_cdf_vec)]), max([0, 1.05 * max(max_cdf_vec)])]);

% Inverse CDF
figure(first_fignum + 2); hold on
for i = 1:N
    plot(xvec, -cdf_diff{i}, 'color', colors(i), 'linestyle', s.line, 'linewidth', s.width);
end
legend(legend_text, 'FontSize', s.fontsize, 'location', 'best'); hold off
min_icdf_vec = zeros(N, 1);
max_icdf_vec = zeros(N, 1);
for i = 1:N
    min_icdf_vec(i) = min(-cdf_diff{i});
    max_icdf_vec(i) = max(-cdf_diff{i});
end
xlim([0, 1]);
xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b')
ylabel('Difference in Inverse CDF', 'FontSize', s.fontsize, 'FontWeight', 'b');
ylim([min([0, 1.05 * min(min_icdf_vec)]), max([0, 1.05 * max(max_icdf_vec)])]);

% Calculate summary statistics
stab100_vec = zeros(N, 1);
stab50_vec = zeros(N, 1);
sum_stat1 = sum_stat_cell{1};
for i = 1:N
    stab100_vec(i) = sum_stat_cell{i + 1}.ergodic.stab_frac - sum_stat1.ergodic.stab_frac;
    stab50_vec(i) = sum_stat_cell{i + 1}.ergodic.stab_frac_pct50 - sum_stat1.ergodic.stab_frac_pct50;
end

% Save plots if desired
if nargin > 6
    if save_plots
	    f = figure(first_fignum);
	    saveas(f, [figurespath, '/pdf_diff.', s.graphics_format]);
	    f = figure(first_fignum + 1);
	    saveas(f, [figurespath, '/cdf_diff.', s.graphics_format]);
	    f = figure(first_fignum + 2);
	    saveas(f, [figurespath, '/inverse_cdf_diff.', s.graphics_format]);
    end
end

% Return output struct
stab_summary_stats.stab100_vec = stab100_vec;
stab_summary_stats.stab50_vec = stab50_vec;

end
