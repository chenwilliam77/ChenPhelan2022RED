function solution = get_calibration_results(etaout, fout, s)
% This function calculates the minimal number of quantities
% for the calibration exercise.
%
% Written by William Chen and Greg Phelan, Jan. 2022
    % Set up
    N                 = length(etaout);
    Dyn               = zeros(N, 8);
    fp                = zeros(N, 3);
    Sharpe            = zeros(N, 1);
    risk_premium      = zeros(N, 1);
    liquidity_premium = zeros(N, 1);

    % Calculate relevant equilibrium quantities quantities
    % Note that dyn in fnct is [psi, sigma_eta, sigma_Q, mu_eta, mu_theta, sigma_theta];
    for n = 1:N
        [fp(n, :), Dyn(n, :), risk_premium(n), Sharpe(n), liquidity_premium(n)] = fnct(etaout(n), fout(n, :), s);
    end

    % Calculate derivatives via FD to get muQ
    Q = fout(:, 3);
    Qp = fp(:, 3);
    Qp(1:end - 1) = diff(Q) ./ diff(etaout); % forward finite difference w/zero boundary condition at end
    Qpp = second_deriv(etaout, Q);
    Qpp(1) = 0; % sigma_eta = 0 at eta = 0
    muQ = Qp ./ Q .* Dyn(:, 4) .* etaout + Qpp ./ Q ./ 2 .* (Dyn(:, 2) .* etaout).^2;

    % Calculate the risk free rate
    % by finding interpolation nodes associated with
    % reasonably valued risk-free rates and using a spline
    iota       = investment_fnct(Q, 'iota', s);
    muK        = investment_fnct(Q, 'Phi', s) - s.delta;
    invst_elas = investment_fnct(Q, 'invst_elasticity', s);
    DS         = deposit_spread(Dyn(:, 1), Dyn(:, 6), etaout, s);
    ab         = arrayfun(@(x) productivity_fnct(x, 0, s), etaout);
    if isfield(s, 'leverage_constraint') % check for whether there are leverage constraints
        drf = zeros(N, 1);
        for i = 1:N
            L = leverage_constraint(etaout(i), fout(i, :), s.leverage_constraint);
            if Dyn(i, 1) >= etaout(i) * (L + 1)
                % leverage constraint binds => use HH FOC: E[drh] - drf = sigma_Ch * ssQ
                % Note sigma_Ch = Dyn(:, 5); sigma_Q = Dyn(:, 3)
                drf(i) = (s.ah - iota(i)) ./ Q(i) + muK(i) + muQ(i) + ...
                    s.sigma * Dyn(i, 3) - Dyn(i, 5) * (s.sigma + Dyn(i, 3));
            else
                % leverage constraint does not bind => use bank FOC
                drf(i) = (ab(i) - iota(i)) ./ Q(i) + muK(i) + muQ(i) + ...
                    s.sigma * Dyn(i, 3) + DS(i) - (risk_premium(i) + liquidity_premium(i));
            end
        end
    else % banks' FOC always holds
        drf   = (ab - iota) ./ Q + muK + muQ + s.sigma * Dyn(:, 3) + DS - (risk_premium + liquidity_premium);
    end
    inds1     = find(abs(drf) < 0.1); % bigger than 0.1 => spike near non-convexity
    inds2     = find(abs(etaout - s.etaPUT) > .005); % typically also spikes in eta nodes
    inds3     = find(abs(etaout - s.etaLAW) > .005); % near etaPUT and etaLAW
    goodinds  = intersect(intersect(inds1, inds2), inds3); % good nodes are those away from these points
    goodinds  = goodinds(1:3:length(goodinds)); % Thin out the interpolation further
    drf       = interp1(etaout(goodinds), drf(goodinds), etaout, 'pchip', 'extrap');

    % Calculate derivatives via FD to get GDP growth rate and volatility
    gdp = (ab - s.ah) .* Dyn(:, 1) + s.ah;
    gdp_p = zeros(length(etaout), 1);
    gdp_p(1:end - 1) = diff(gdp) ./ diff(etaout);
    gdp_pp = second_deriv(etaout, gdp);
    gdp_pp(1) = 0; % sigma_eta = 0 at eta = 0
    mu_gdp = gdp_p ./ gdp .* Dyn(:, 4) .* etaout + gdp_pp ./ gdp ./ 2 .* (Dyn(:, 2) .* etaout).^2 + ...
        muK + (gdp_p ./ gdp .* Dyn(:, 2) .* etaout) .* s.sigma;
    sigma_gdp = (gdp_p ./ gdp .* Dyn(:, 2) .* etaout) + s.sigma;

    % Calculate derivatives via FD to get investment growth rate and volatility
    iota_p = zeros(length(etaout), 1);
    iota_p(1:end - 1) = diff(iota) ./ diff(etaout);
    iota_pp = second_deriv(etaout, iota);
    iota_pp(1) = 0; % sigma_eta = 0 at eta = 0
    mu_iota = iota_p ./ iota .* Dyn(:, 4) .* etaout + iota_pp ./ iota ./ 2 .* (Dyn(:, 2) .* etaout).^2 + ...
        muK + (iota_p ./ iota .* Dyn(:, 2) .* etaout) .* s.sigma;
    sigma_iota = (iota_p ./ iota .* Dyn(:, 2) .* etaout) + s.sigma;

    % Calculate interest rates
    N = length(etaout);
    interestvec = zeros(N, 1);
    for n = 1:N
        interestvec(n) = interest(etaout(n), s);
    end

    % Calculate ergodic density
    [eta_density, density, cdf, ~] = ergocalc(etaout, Dyn, s);

    % Get inflation
    inflation = interestvec - drf;

    % Populate solution struct
    solution.eta = etaout;
    solution.psi = Dyn(:, 1);
    solution.MarketLeverageRatio = (Dyn(:, 1) ./ etaout - 1) ./ fout(:, 1) + 1;
    solution.sigma_eta = Dyn(:, 2);
    solution.mu_eta = Dyn(:, 4);
    solution.sigma_Ch = Dyn(:, 5);
    solution.Sharpe = Sharpe;
    solution.drf = drf;
    solution.deposit_spread = DS;
    solution.iota = iota;
    solution.gdp = gdp;
    solution.sigma_gdp = sigma_gdp;
    solution.interestvec = interestvec;
    solution.inflation = inflation;
    solution.eta_density = eta_density;
    solution.density = density;
    solution.invst_elasticity = invst_elas;
end
