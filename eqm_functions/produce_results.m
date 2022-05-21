function [solution, summary_stats] = produce_results(etaout, fout, s)
% [solution, summary_stats] = produce_results(etaout, fout, s)
%
% This function calculates equilbrium variables, computes summary statistics,
% and plots results, given solutions for q, and
% theta of the ODE.  solution.eta is grid for eta; fout is matrix for
% theta, theta', Q. s is structure of paramaters.
%
% Written by William Chen and Gregory Phelan, Jan. 2022

% Calculate results
if isfield(s, 'calibrate') & s.calibrate
    % Compute only necessary moments for calibration

    solution = get_calibration_results(etaout, fout, s);

    % Calculate summary statistics
    summary_stats = get_calibration_summary_stats(solution, s);
else
    solution = get_results(etaout, fout, s);

    % Calculate summary statistics
    summary_stats = get_summary_stats(solution, s);

    if isfield(s, 'plot_results') & s.plot_results
        plot_results(solution, s);
    end
end
