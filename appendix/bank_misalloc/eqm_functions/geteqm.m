% This function solves for equilibrium. Takes one input, where `s` is a struct
% of settings (e.g. parameters, numerical settings).
%
% written by  William Chen and Greg Phelan, Jan. 2021
function [solution, summary_stats] = geteqm(s)

    if isfield(s, 'verbose')
        if s.verbose
            tic
        end
    end

    %% Set up

    % Calculate some numbers for later computations
    s.Qbeg = consumption_market_clearing(0, 1, 0, 'Q', s);
    s.Qend = consumption_market_clearing(1, 1, 1, 'Q', s);

    % initial conditions for (theta, thetap, Q) at eta = 0 or eta = eta_TR
    % Search over theta and/or thetap to
    % solve for equilibrium. If search_init_conds = 1,
    % then we search over theta,
    % and if search_init_conds = 2,
    % then we search over both theta and thetap.
    % In both cases, we start with a bisection search
    % over theta, and then we solve the nonlinear root
    % problem using a minimization algorithm.
    if isfield(s, 'search_init_conds')
        search_init_conds = s.search_init_conds;
    else % Default to searching over 1 initial condition
        search_init_conds = 1;
    end
    theta0  = 1;                     % just initializing
    thetap0 = s.thetap0;             % A large negative number
    Q0      = s.Qbeg * s.Q0_perturb; % Q(0) = (ah + deprec_adjust / eps2) / (r + eps1 / eps2)

    etaspan = [s.start s.end]; % we start just off zero because of the singularity there
    F0 = [theta0; thetap0; Q0];

    % Set up termination events, see `evntfcn.m` for details
    odeevnt = @(eta, f) evntfcn(eta, f, s);

    % Set up ODE solver
    options = odeset('events', odeevnt, 'AbsTol', s.abs_tol, 'RelTol', s.rel_tol);
    odefun = @(eta, f) fnct(eta, f, s);
    if isfield(s, 'autodiff_jacobian') % use automatic differentiation for ODE
        if s.autodiff_jacobian         % (but since ODE has 3 variables, unlikely to speed up)
            Jf = @(eta, f, yp) full(getderivs(odefun(eta, myAD(f))));
            options = odeset(options, 'Jacobian', Jf);
        end
    end
    if isfield(s, 'odesolver') % Choose solver
        if strcmp(s.odesolver, 'ode45')
            odesolver = @(F0) ode45(odefun, etaspan, F0, options);
        elseif strcmp(s.odesolver, 'ode15s')
            odesolver = @(F0) ode15s(odefun, etaspan, F0, options);
        elseif strcmp(s.odesolver, 'ode23')
            odesolver = @(F0) ode23(odefun, etaspan, F0, options);
        elseif strcmp(s.odesolver, 'ode23s')
            odesolver = @(F0) ode23s(odefun, etaspan, F0, options);
        else
            error('Only ode45, ode23, ode23s, and ode15s are currently permitted.');
        end
    else % Default to ode45
        odesolver = @(F0) ode45(odefun, etaspan, F0, options);
    end

    % the true ODE is a boundary-value problem, which is much easier to solve as
    % an initial-value problem (IVP) and iterate until boundary conditions
    % are matched. See the paper's appendix for details on the IVP being solved.
    thetaL = s.thetastar; thetaR = s.theta0R;
    for iter = 1:s.bisect_iter % Begin with bisection to narrow where theta should be
        F0(1) = (thetaL + thetaR) / 2;

        [~, ~, ~, ~, IE] = odesolver(F0);

        if any(IE == 1)         % theta hit 1 too early
            thetaL = F0(1);
        elseif ~isempty(IE)     % thetap or Qp hit zero but theta > 1 still
            thetaR = F0(1);
        else
            keyboard; % initiate debugging
        end
    end

    if isfield(s, 'lsq_iters')
        optimoptions = optimset('MaxIter', s.lsq_iters, 'display', 'off');
    else
        optimoptions = optimset('display', 'off');
    end

    % Figure out which solver to use
    if strcmp(s.lsq_solver, 'none')
        % Nothing happens, we use (thetaL + thetaR) / 2 as the final
        % guess for theta0
    elseif strcmp(s.lsq_solver, 'lsqnonlin')
        if isfield(s, 'lsq_boundaries')
            a = s.lsq_boundaries(1);
            b = s.lsq_boundaries(1);
        else
            a = thetaL;
            b = thetaR;
        end

        [out, resnorm, ~, exitflag] = lsqnonlin(@(x) fit_boundary_conditions(x, F0, s), (a + b) / 2, a, b, optimoptions);
    elseif search_init_conds == 1
        if strcmp(s.lsq_solver, 'fminsearch')
            [out, resnorm, exitflag] = fminsearch(@(x) fit_boundary_conditions(x, F0, s), (thetaL + thetaR) / 2, optimoptions);
        elseif strcmp(s.lsq_solver, 'fmincon')
            [out, resnorm, exitflag] = fmincon(@(x) fit_boundary_conditions(x, F0, s), (thetaL + thetaR) / 2, ...
                [], [], [], [], thetaL, thetaR, [], optimoptions);
        else
           error(['Cannot use solver ', s.lsq_solver]);
        end
    elseif search_init_conds == 2
        if strcmp(s.lsq_solver, 'fminsearch')
            [out, resnorm, exitflag] = fminsearch(@(x) fit_boundary_conditions(x, F0, s), ...
                [(thetaL + thetaR) / 2; F0(2)], optimoptions);
        elseif strcmp(s.lsq_solver, 'fmincon')
            if isfield(s, 'thetap0_bounds')
                thetap0_min = min(s.thetap0_bounds);
                thetap0_max = max(s.thetap0_bounds);
                [out, resnorm, exitflag] = fmincon(@(x) fit_boundary_conditions(x, F0, s), ...
                    [(thetaL + thetaR) / 2; F0(2)], [], [], [], [], ...
                    [thetaL; thetap0_min], [thetaR; thetap0_max], [], optimoptions);
            else
                error(['Cannot use solver fmincon without specifying thetap0_bounds in settings struct']);
            end
        else
           error(['Cannot use solver ', s.lsq_solver]);
        end
    else
        error(['search_init_conds in settings struct must be 1 or 2']);
    end

    % Calculate final solution
    if search_init_conds == 1
        F0(1) = out;
    else
        F0(1) = out(1);
        F0(2) = out(2);
    end
    [etaout, fout, ~, ~, ~] = odesolver(F0);

    % Check if solution satisfies error tolerances
    if resnorm > s.lsq_tol
        if s.verbose
            disp(fout(end, :));
            disp(['Exit flag for the minimization algorithm is ', num2str(exitflag)]);
        end
        error(['Unable to find a minimum with a residual norm below ', num2str(s.lsq_tol), ...
            '. The current residual norm is ', num2str(resnorm)]);
    end

    % Produce results and plot them.
    [solution, summary_stats] = produce_results(etaout, fout, s);

    % Print out summary statistics of interest
    if isfield(s, 'verbose')
        if s.verbose
            if isfield(s, 'crisis_output_losses') && s.crisis_output_losses
                avg_gdp_loss = crisis_output_losses(solution, 5, 1000, s);
            end
	        etaPUT_psi = interp1(solution.eta, solution.psi, s.etaPUT, s.interp);
            if s.etaLAW < solution.eta(end)
                etaLAW_psi = interp1(solution.eta, solution.psi, s.etaLAW, s.interp);
            else
                etaLAW_psi = 1; % large number b/c etaLAW is too large in this case
            end

            toc;
            fprintf('\n');

            fprintf('The ergodic means for calibration targets are ...\n\n');
            fprintf('Prob. of Crisis:            %.2f%%\n', 100 - summary_stats.ergodic.stab_frac_pct50);
            fprintf('GDP Vol:                    %.2f%%\n', 100 * summary_stats.ergodic.sigma_gdp);
            fprintf('Risk-Free Rate:             %.2f%%\n', 100 * summary_stats.ergodic.drf);
            fprintf('Market Leverage Ratio:      %.2f\n', summary_stats.ergodic.MarketLeverageRatio);
            fprintf('Normal Times psi:           %.2f%%\n', 100 * etaLAW_psi);
            fprintf('Crisis psi:                 %.2f%%\n', 100 * etaPUT_psi);
            fprintf('Inflation not in ZLB:       %.2f%%\n', 100 * summary_stats.ergodic.non_zlb_pi);

            fprintf('\n');
            fprintf('The ergodic means for some untargeted moments are ... \n\n')
            fprintf('Sharpe:                     %.2f%%\n', 100 * summary_stats.ergodic.Sharpe);
            fprintf('Investment Price Elas.:     %.2f\n', summary_stats.ergodic.invst_elasticity);
            if isfield(s, 'crisis_output_losses') && s.crisis_output_losses
                fprintf('Output Loss in Crisis:      %.2f%%\n', 100 * avg_gdp_loss);
            end
            fprintf('Investment-Capital Ratio:   %.2f%%\n', 100 * summary_stats.ergodic.iota);
            fprintf('Prob. of Distress:          %.2f%%\n', 100 - summary_stats.ergodic.stab_frac);
        end
    end
end
