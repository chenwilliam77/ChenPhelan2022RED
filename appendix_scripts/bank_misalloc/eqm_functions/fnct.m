function [fp, dyn, risk_premium, Sharpe, liquidity_premium] = fnct(eta, f, s)
% takes scalar eta, a 3x1 vector f = [theta, theta', Q], and parameters
% s, and computes fp, the derivative of f with respect to eta.
% Capital growth is driven by internal investment.
%
% written by William Chen and Greg Phelan, Jan. 2022

	% Set up
	theta  = f(1);
	thetap = f(2);
	Q      = f(3);

	if s.nomrule == 0
	    int = 0;
	else
	    int = interest(eta, s);
	end
	% liquidity_premium = s.kappa / s.m * (int - s.itarget); % If using this, will need to check that other uses in this script of liquidity_premium still hold
    liquidity_premium = s.kappa / s.m * int;
    ab = productivity_fnct(eta, 0, s);

	% Market clearing for consumption pins down psi
    psi = consumption_market_clearing(eta, theta, Q, 'psi', s);
	% psi = (Q * (s.r * (1 + (theta - 1) * eta) + s.eps1 / s.eps2) - ...
    %     s.ah - (exp(s.eps2 * s.delta) - s.eps2 * s.delta) / s.eps2) / (s.ab - s.ah);

    % Calculate share of household wealth held in non-bank net worth
	% omega = nh / wh = nh / (nh + theta nb)
    %       = (nh / (QK)) / ((nh + theta nb) / (QK))
    %       = (1 - eta) / ((1 - eta) + theta eta)
	omega = (1 - eta) / ((1 - eta) + theta * eta);

    % Compute deposit spread, given psi, omega, and eta
    DS = deposit_spread(psi, omega, eta, s);

    % Check if leverage constrained
    if isfield(s, 'leverage_constraint')
        [L, dL] = leverage_constraint(eta, f, s.leverage_constraint);
    else
        L = Inf;
    end

    if psi >= eta * (L + 1)
        psi = eta * (L + 1);
        lc = true;
    else
        lc = false;
    end

    % Calculate ssq and Qp
    if psi >= 1 % Enforce that psi <= 1, and calculate Qp analytically
        psi     = 1;
        ab_p    = productivity_fnct(eta, 1, s);
        Qp      = (ab_p - Q * s.r * (thetap * eta + (theta - 1))) / ...
            (s.r * (1 + (theta - 1) * eta) + s.eps1 / s.eps2);
        sigma_Q = Qp * (psi - eta) * s.sigma / (Q - Qp * (psi - eta));
        ssQ     = s.sigma + sigma_Q;
    elseif lc
        % Obtain Qp directly from market-clearing for consumption
        ab_p    = productivity_fnct(eta, 1, s);
        Qp      = (ab_p * eta * (L + 1) + (ab - s.ah) * (L + 1 + eta * dL) - ...
            Q * s.r * (thetap * eta + theta - 1)) / ...
            (s.r * (1 + (theta - 1) * eta) + s.eps1 / s.eps2);
        sigma_Q = Qp * (psi - eta) * s.sigma / (Q - Qp * (psi - eta));
        ssQ     = s.sigma + sigma_Q;
    else % Apply asset pricing condition to solve for ssq = sigma + sigma_q since psi < 1
        ssQ     = sqrt(((ab - s.ah) / Q - liquidity_premium + DS) * ...
              theta / -thetap / (psi - eta));
        sigma_Q = ssQ - s.sigma;
        Qp      = Q * sigma_Q / ((psi - eta) * ssQ);
    end

    % Calculate leverage
    xh    = (1 - psi) / (1 - eta);                 % household leverage
    xb    = psi / eta;                             % bank leverage

    % Calculate volatilities of eta, theta, and household consumption
    sigma_eta     = (psi / eta - 1) * ssQ;
	eta_sigma_eta = sigma_eta * eta;
    sigma_theta   = eta_sigma_eta * thetap / theta;
    sigma_ch      = (omega * (xh - xb) + xb) * ssQ + (1 - omega) * sigma_theta;

    % Calculate drifts
	T           = s.T;
    transfers_b = -T + (psi - eta) * liquidity_premium;
    transfers_h = (eta * T / (1 - eta)) + (psi - eta) * liquidity_premium; % transfers to household per dollar of non-bank equity
    iota        = investment_fnct(Q, 'iota', s);
    if lc
        mu_theta = transfers_h - transfers_b - (xb - 1) * DS - liquidity_premium + sigma_ch * sigma_theta - ...
            xb * ((ab - s.ah) / Q - (-sigma_theta * ssQ + liquidity_premium));
        mu_eta   = xb * ab / Q - (xb - 1) * s.ah / Q - iota / Q ...
            - (xb - 1) * liquidity_premium + (xb - 1) * DS + transfers_b ...
            + (xb - 1) * (sigma_ch - ssQ) * ssQ;
    else
        mu_theta = transfers_h - transfers_b - liquidity_premium + DS + sigma_ch * sigma_theta;
        mu_eta   = (ab - iota) / Q + transfers_b + (xb - 1) * ((sigma_ch - sigma_theta) - ssQ) * ssQ;
    end

    % Calculate differential equation for thetapp
    eta_mu_eta = eta * mu_eta;
    thetapp    = 2 * (mu_theta * theta - thetap * eta_mu_eta) / eta_sigma_eta^2;

    fp = [thetap; thetapp; Qp];

    % Calculate remaining values of interest
    if lc
        % excess_returns = E[dr_b] - dr_d = (ab - s.ah) / Q + sigma_ch * ssQ + DS
        % and E[dr_b] - dr_d = risk_premium + liquidity_premium
        risk_premium = (ab - s.ah) / Q + sigma_ch * ssQ + DS - liquidity_premium;
    else % HH risk premium + excess risk premium
        risk_premium = (sigma_ch - sigma_theta) * ssQ;
    end
    Sharpe           = (risk_premium + liquidity_premium) / ssQ; % (E[dr_b] - dr_d) / vol[dr_b]

    dyn = [psi, sigma_eta, sigma_Q, mu_eta, sigma_ch, omega, mu_theta, sigma_theta];
end
