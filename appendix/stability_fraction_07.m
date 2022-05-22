% This script plots the Probability of Distress
% under a discontinuous piecewise rule with strike at 3%
%
% Written by William Chen and Gregory Phelan, Mar. 2022

clear;
close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

figurespath = [pwd(), '/../figures/stab_frac_07'];
if ~isfolder(figurespath)
    mkdir(figurespath);
end
save_plots = 1; % set to 1 to save plots

calibration_parameters;
s.plot_results = 0;
s.verbose = 0;
s.nomruletype = 1; % switch to piecewise
set(groot, 'DefaultAxesFontSize',  s.fontsize); % Change default axes fonts
set(groot, 'DefaultTextFontSize',  s.fontsize); % Change default text fonts
set(groot, 'DefaultLineLineWidth', s.width);

colors = get(gca, 'ColorOrder');
iter = 0;

%% Get Stability fraction, varying piecewise int rules
% using piecewise rules with i^LAW (high eta) and i^Put (low eta),
% calculating stability on the grid of policy rules

% Set monetary policy
N1 = 3;
tic
s.center = .07;
s.etaPUT = s.center;
s.etaLAW = s.center;

i1max  = .03;   % PUT max
i2max  = .05;   % LAW max
s.imax = i2max;

int1grid = linspace(0, i1max, N1);
% int2grid = 0.0:.01:i2max;
int2grid = 0:(5/3/100):0.05;
% int2grid = [0, 0.02, 0.03, 0.04, 0.05];
N2 = length(int2grid);

mat1 = repmat(int1grid', 1,  N2);
mat2 = repmat(int2grid,  N1, 1);
ss = cell(N1 * N2, 1);
v1 = reshape(mat1, N1 * N2, 1);
v2 = reshape(mat2, N1 * N2, 1);
v0 = zeros(N1 * N2, 1);

for i = 1:(N1 * N2)
    ss{i} = s; % Copy parameters struct so no problems with changing s.
end
err_msg = cell(N1, N2);
% parfor n = 1:(N1 * N2)
for n = 1:(N1 * N2)
    try
        v0(n) = get_stab_frac([v1(n); v2(n)], ss{n}, 'mp');
    catch err
        err_msg{n} = err;
    end
end
Vgrid = reshape(v0, [N1, N2]);

toc

% Count number of errors
err_ct = 0;
for i = 1:numel(err_msg)
    if ~isempty(err_msg{i})
        err_ct = err_ct + 1;
    end
end
if err_ct == 0
    fprintf('No errors encountered during calculation of stability fractions.\n');
else
    fprintf(['Errors occurred for ', num2str(err_ct), ' policy rules.\n']);
end

figure(41);
plot(100 * int2grid, 100 - Vgrid);

legstr = '$i^{Put}=';
legendInfo = cell(N1, 1);
for n = 1:N1
    legendInfo{n} = strcat(legstr, num2str(100 * int1grid(n), 2), '$\%');
end
% legend(legendInfo, 'location', 'best') % use legend in stability_fraction_03.m

xlabel('$i^{LAW}$\%')
ylabel('Prob. of Distress')

if save_plots
    f = figure(41);
    saveas(f, [figurespath, '/stability_v_ilaw.', s.graphics_format]);
end

close 1
