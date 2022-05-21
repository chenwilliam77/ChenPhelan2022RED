% This script maximizes expected welfare for a Fed Put.
% Thus, there are 4 arguments over which to
% optimize, iput, ilaw, etaPUT, and etaLAW - etaPUT.m
%
% Written by William Chen and Gregory Phelan, Feb. 2022

close all;
addpath ../parameters_specs
addpath(genpath('../eqm_functions'))

savepath = [pwd(), '/../save/output_data/max_expV_fedput_4args/'];
if ~isdir(savepath)
    mkdir(savepath);
end
run_fmincon = 1; % set to 1 to minimize the function over the hyper cube specified by lb and ub
nworkers    = 16;

% Set up
calibration_parameters;
s.plot_results = 0;
s.verbose = 0;
s.max_expV = 1;
thetapgrid = [s.thetap0, -1e3, -5e3, -1e4, -1e5, -5e5, -1e6, -1e8];

% fmincon settings
lb = [-eps; -eps; 0.0; 0.0]; % eps to allow zero as possible value
ub = [.04; s.imax; .05; .10];
algo = 'interior-point';
Ns = [5, 5, 5, 5]; % dimension of initial guesses w/in hyper cube specified by lb and ub
lb_pad = 1e-5;  % padding to ensure guesses are w/in interior of hyper cube
ub_pad = -1e-5;

% Create grid of initial guesses
inputs = cell(Ns);
ss = cell(Ns);
x0_1D = cell(length(Ns));
RHS_str = '';
LHS_str = '';
input_str = cell(length(Ns));
for i = 1:length(Ns)
	x0_1D{i} = linspace(lb(i) + lb_pad, ub(i) - ub_pad, Ns(i))';
	RHS_str   = [RHS_str, 'x0_1D{', num2str(i), '},'];
	LHS_str   = [LHS_str, 'x0_ndgrid', num2str(i), ','];
	input_str{i} = ['x0_ndgrid', num2str(i)'];
end
RHS_str = RHS_str(1:end - 1); % remove final comma
LHS_str = LHS_str(1:end - 1); % remove final comma
eval(['[', LHS_str, ']=ndgrid(', RHS_str, ');']); % create ndgrids for each element of input guess
for i = 1:numel(inputs) % populate inputs
    ss{i} = s; % Copy struct s
    joint_str = '';
    for j = 1:length(Ns)
        joint_str = [joint_str, 'x0_ndgrid', num2str(j), '(', num2str(i), ');'];
    end
    joint_str = joint_str(1:end - 1); % remove final semi-colon
    eval(['inputs{i} = [', joint_str, '];']);
end

% Output arguments
best_rules = cell(Ns);
fvals = zeros(Ns);
exitflags = cell(Ns);

% Run fmincon!
if run_fmincon
    options = optimoptions('fmincon', 'Algorithm', algo);
    myprocs = parpool(nworkers);
    parfor i = 1:numel(inputs)
        x0 = inputs{i};
	    [best_rules{i}, fvals(i), exitflags{i}, ~] = fmincon(@(x) get_ExpV_maxfnct(x, ss{i}, 'fp', 'thetap0', thetapgrid, 'nan_error', 1, 'min_problem', 1), ...
	        x0, [], [], [], [], lb, ub, [], options);
    end
	delete(myprocs);

    [best_fval, best_rule_i] = min(fvals, [], 1:length(Ns), 'linear'); % min b/c minimization problem
	best_rule = best_rules{best_rule_i};

    save([savepath, '/fmincon_run_', algo, '_', date, '.mat'], 'lb', 'ub', 'inputs', 'ss', 'thetapgrid', 'best_rules', 'fvals', 'exitflags', 'best_rule', 'best_fval');
end

exit
