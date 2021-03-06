function out = fit_boundary_conditions(x, F0, s)
% out = fit_boundary_conditions(x, F0, s)
%
% Wrapper function reformulating the
% boundary value problem as a least squares problem
% by minimizing the residuals of the boundary conditions
% at etastar.
%
% Written by William Chen and Gregory Phelan, Jan. 2022

    % Determine left boundary conditions
    if isfield(s, 'tail_risk')
        tail_risk = s.tail_risk;
        if tail_risk && ~isfield(s, 'eta_TR')
            error('Cannot apply tail risk insurance if eta_TR is not set.');
        end
    else
        tail_risk = false;
    end

    % Set up termination events, see `evntfcn.m` for details
    odeevnt = @(eta, f) evntfcn(eta, f, s);

    % Set up ODE solver
    if length(x) == 1 % if input has length 1, then we are searching over theta
        F0(1) = x;
    elseif length(x) == 2 % if input has length 2, then we are searching over thetap
        F0(1) = x(1);
        F0(2) = x(2);
    else
        error('Length of input x to fit_boundary_conditions must be 1 or 2.');
    end
    options = odeset('events', odeevnt, 'AbsTol', s.abs_tol, 'RelTol', s.rel_tol);
    odefun = @(eta, f) fnct(eta, f, s);
    if isfield(s, 'autodiff_jacobian')
        if s.autodiff_jacobian
            Jf = @(eta, f, yp) full(getderivs(odefun(eta, myAD(f))));
            options = odeset(options, 'Jacobian', Jf);
        end
    end
    if tail_risk
        % Calculate initial conditions
        F0(2) = -F0(1) / s.eta_TR; % implication of tail risk insurance
        tr_options = optimset('display', 'off');
        tr_out = fsolve(@(y) tail_risk_init_cond(y, F0(1), s), [(s.Qmax + s.Qmin) / 2;  s.eta_TR + 1e-3], tr_options);
        F0(3) = tr_out(1); % implied Q satisfying Q' = 0.

        % Set up system of equations for initial condition
        etaspan = [s.eta_TR s.end];
    else
        etaspan = [s.start s.end]; % we start just off zero because of the singularity there
    end

    if isfield(s, 'odesolver')
        if strcmp(s.odesolver, 'ode45')
            odesolver = @(F0) ode45(odefun, etaspan, F0, options);
        elseif strcmp(s.odesolver, 'ode15s')
            odesolver = @(F0) ode15s(odefun, etaspan, F0, options);
        elseif strcmp(s.odesolver, 'ode23')
            odesolver = @(F0) ode23(odefun, etaspan, F0, options);
        elseif strcmp(s.odesolver, 'ode23s')
            odesolver = @(F0) ode23s(odefun, etaspan, F0, options);
        else
            error('Only ode45, ode23, ode23s, and ode15s are currently permitted as solvers.');
        end
    else % Default to ode45
        odesolver = @(F0) ode45(odefun, etaspan, F0, options);
    end

    % Solve ODE
    [etaout, fout, ~, ~, ~] = odesolver(F0);

    % Construct residuals
    [value, ~, ~] = odeevnt(etaout(end), fout(end, :));
    residuals = [value'; (fout(end, 3) - s.Qmax)];

    if strcmp(s.lsq_solver, 'lsqnonlin')
        out = s.lsq_weights .* residuals;
    elseif strcmp(s.lsq_solver, 'fminsearch') || strcmp(s.lsq_solver, 'fmincon')
        out = sum(s.lsq_weights .* residuals.^2);
    else
        error(['Cannot use solver ', s.lsq_solver]);
    end
end
