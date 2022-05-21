% This script defines the parameters and numerical settings
% to solve the model under Calibration

% Preferences
s.r         = 0.02;     % time preference
s.lambda    = 419.4175; % inflation costs
s.pi_target = 0.02;     % target inflation rate

% Convenience yield parameters
% v(x_hd) = beta_1 * (x_hd + beta_2) ^ (1 - beta_3) / (1 - beta_3)
s.beta_1 = 0;
s.beta_2 = 0;
s.beta_3 = 0;

% Technology
s.ab    = 0.127;       % bank productivity
s.ah    = 0.8 * s.ab;  % HH productivity
s.eps1  = 1.0;         % capital efficiency parameter (see investment_fnct.m for details)
s.eps2  = 3.0;         % capital adjustment costs (this is the epsilon reported in the paper, see investment_fnct.m for details)
s.T     = .048;        % tax on banks' wealth
s.delta = .10;         % Depreciation rate
s.sigma = 0.015;       % volatility of capital growth shocks

% Monetary Policy parameters
s.chi   = 0.10;  % frequency of aggregate deposit withdrawal shocks
s.rho   = 0.21;  % fire sale discount
s.m     = 4;     % reserve multiplier, from Drechsler et al. (2018)
s.kappa = .4085; % exposure of deposits to fire sales, from Drechsler et al. (2018)
s.imax  = s.chi * s.rho / (1 - s.rho) * s.kappa / (1 + s.kappa) / (s.kappa / s.m);

% Type of monetary policy rule
s.nomrule     = 1;  % set to 1 to turn on MP rule, set to 0 to have i = 0 everywhere
s.nomruletype = 0;  % set to 1 for piecewise, else polynomial with bounds
s.massdiff    = 0;

% Piecewise MP rule, two break points
if s.nomruletype
    s.etaLAW = .1;
    s.etaPUT = .05;
    s.ilaw = 0.02;
    s.ireg = 0.02;
    s.iput = 0;
else
    % Linear MP rule, two break points
    s.etaPUT = .02; % start cutting once crisis region is hit.
    s.etaLAW = s.etaPUT + .04; % start cutting once financial distress appears
    s.ilaw = 0.041; % Target 2% inflation outside of the zlb, also modal steady state nominal interest rate is 4% in FRBNY DSGE
    s.ireg = 0; % in case you want to switch policy rules
    s.iput = 0; % ZLB
end

% Boundary conditions (affects boundary macropru)
s.issuance = Inf; % = 1 + gamma, where gamma is the cost of issuing new equity
s.thetastar = 1;  % Boundary condition for theta at etastar

% for plots
s.color = 'k'; s.line = '-'; s.width = 3;
s.fontsize = 20; s.use_title = 1;
s.legend_fontsize = 16; s.use_title = 1;
s.graphics_format = 'png';
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', s.fontsize)
set(groot, 'defaultLegendFontSize', s.legend_fontsize)
set(groot, 'defaultLineLineWidth', s.width);

% Numerical parameters
s.start = 1e-3; % left boundary of eta grid
s.end = 1 - 1e-3; % maximum right boundary of eta grid
s.interp = 'linear'; % interpolation method for geteqm (e.g. integration for ergodic moments)
s.thetap0 = -5e4; % Smaller guesses (in magnitude) tend to fix problem of non-monotonic psi, and not matching boundary conditions
s.theta0L = 1; % lowest theta0 for initial bisection search
s.theta0R = 1e3; % largest theta0 for initial bisection search
s.Q0_perturb = 1.000; % perturbation of Q away from its theoretical minimum at eta = s.start
s.rel_tol = 1e-6; % relative tol for ODE
s.abs_tol = 1e-3 * s.rel_tol; % absolute tol for ODE (MATLAB default sets abs_tol = 1e-3 * rel_tol)
s.bisect_iter = ceil((s.theta0R - s.theta0L) / 200); % number of bisection iterations, e.g. if s.theta0R = 1e3, s.theta0L = 1, then this equals 5
s.verbose = 1; % print output of interest (e.g. was convergence achieved? summary statistics?)
s.plot_results = 1;
s.N_welfare = 100; % grid points for calculating welfare via finite differences
s.odesolver = 'ode45';
s.lsq_solver = 'fminsearch'; % preferred for now, seems to work better b/c lsq problem is not smooth, and fminsearch is derivative free
s.lsq_weights = [20; 1; 20]; % weights on boundary conditions [theta(etastar); theta'(etastar); Q(etastar)];
s.lsq_iters = 100; % maximum number of iterations for lsq solver
s.lsq_tol = 1e-2; % tolerance level for lsq solver, 1e-2 is sufficient to deliver good results for boundary condition
s.calibrate = 0; % if set to 1, then the solvers are run in "calibration" mode, with fewer quantities being calculated to accelerate computation
s.Phi_form = 1; % functional form for investment, see investment_fnct.m
