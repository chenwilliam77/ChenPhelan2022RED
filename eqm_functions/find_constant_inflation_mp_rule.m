function [new_s, converged] = find_constant_inflation_mp_rule(s, varargin)
% [new_s, converged] = find_constant_inflation_mp_rule(s, varargin)
%
% Given other parameters, find a monetary policy rule
% which delivers an approximately constant inflation rate
% specified by `s.pi_target`. Returns a new settings struct
% with the new policy rule.
%
% The policy rule is updated in a tatonnement style.
% At a series of collocation nodes, we aim to set
% inflation close to the target and interpolate
% between these collocation nodes. If the error is nonzero,
% then we update the interest rate at that node toward zero
% with some damping on the magnitude of the step.
%
% There are a number of numerical settings that affects
% how the policy rule is updated and are entered
% as optional arguments via varargin.
%
% interp:        interpolation method
%
% thinning:      an integer which controls
%                how dispersed the nodes are for the constructed
%                policy rule. A higher `thinning` value
%                results in a lower dimensional policy rule.
%
% max_iters:     maximum number of iterations for updating the
%                the monetary policy rule
%
% damping:       damping weight for updating of the policy rule
%
% tol:           convergence tolerance (maximum absolute deviation
%                from target inflation at nodes)
%
% Written by William Chen and Gregory Phelan, Jan. 2022

% Process varargin
for j = 1:2:length(varargin)
    if strcmp(varargin{j}, 'max_iters')
        max_iters = varargin{j + 1};
    elseif strcmp(varargin{j}, 'interp')
        interp = varargin{j + 1};
    elseif strcmp(varargin{j}, 'damping')
        damping = varargin{j + 1};
    elseif strcmp(varargin{j}, 'thinning')
        thinning = varargin{j + 1};
    elseif strcmp(varargin{j}, 'tol')
        tol = varargin{j + 1};
    end
end
if ~exist('max_iters', 'var')
    max_iters = 100; % number of iterations to find constant inflation rate policy rule
end
if ~exist('interp', 'var')
    interp = 'linear'; % linear interpolation
end
if ~exist('damping', 'var')
    damping = 0.5;
end
if ~exist('thinning', 'var')
    thinning = 10; % collocation node every `thinning`th point
end
if ~exist('tol', 'var')
    tol = 1e-2; % convergence tolerance (maximum absolute error)
end

% Operate on a copy of the original settings struct
new_s = s;

% Find initial equilibrium to get a guess of the policy rule
% using as an initial guess a constant nominal rate
% rule equaling the discount rate plus the target inflation rate
new_s.calibrate = 1;
new_s.plot_results = 0;
new_s.verbose = 0;
new_s.nomrule = 1;
new_s.nomruletype = 1;
new_s.iput = s.pi_target + s.r;
new_s.ireg = s.pi_target + s.r;
new_s.ilaw = s.pi_target + s.r;
[sol, ~] = geteqm(new_s);
new_s.nomruletype = s.nomruletype; % set it back to the original Fed Put type

% Update settings struct
new_s.nomrule = 2; % use a user-specified interpolant rather than the interest.m function

% Run update loop
converged = 0;
for iter = 1:max_iters
    [new_s, max_abs_err] = update_policy_rule(sol, new_s, thinning, damping, interp);

    if max_abs_err < tol % break loop if tolerance is met
        if s.verbose
            disp('Convergence achieved. Constant inflation MP rule found.');
        end
        converged = 1;
        break;
    else % Resolve equilibrium otherwise
        if s.verbose
            disp(['Maximum absolute error at iter ' num2str(iter) ': ' num2str(max_abs_err)]);
        end
        [sol, ~] = geteqm(new_s);
    end
end

% Return settings of the struct back to the original
% except for the fact that the MP rule is exogenously
% specified by an interpolant
if isfield(s, 'iput')
    new_s.iput = s.iput;
end
if isfield(s, 'ireg')
    new_s.ireg = s.ireg;
end
if isfield(s, 'ilaw')
    new_s.ilaw = s.ilaw;
end
if isfield(s, 'calibrate')
    new_s.calibrate = s.calibrate;
else
    new_s.calibrate = 0;
end
if isfield(s, 'plot_results')
    new_s.plot_results = s.plot_results;
else
    new_s.plot_results = 0;
end
if isfield(s, 'verbose')
    new_s.verbose = s.verbose;
else
    new_s.verbose = 0;
end
% new_s.nomrule = 2; % keep new_s.snomrule as 2 => interpolant

end

function [new_s, max_abs_err] = update_policy_rule(sol, s, thinning, damping, interp)
    new_s = s;

    % Find the errors based on the thinning parameter
    N = length(sol.eta);
    inds_interp = 1:thinning:N;
    pi_error = s.pi_target - sol.inflation(inds_interp); % deviations in inflation rate from target
    max_abs_err = max(abs(pi_error));

    % Update the policy rule with damping
    i_interp = sol.interestvec(inds_interp) + (1-damping) .* pi_error; % if pi_error < 0, then need to decrease interest rate there
    new_s.nomrule_mp_interp = griddedInterpolant(sol.eta(inds_interp), i_interp, interp, 'linear'); % linear extrapolation
end
