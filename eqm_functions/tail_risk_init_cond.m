function F = tail_risk_init_cond(x, theta, s)
% This function computes residuals to
% pin down the initial condition when
% a tail risk macroprudential policy
% is being used.
%
% Written by William Chen and Greg Phelan, Jan. 2021

Q = x(1);
psi = x(2);
F = reshape(x, 2, 1);

% Set up for residuals
eta = s.eta_TR;
if s.nomrule == 0
    int = 0;
else
    int = interest(eta, s);
end
liquidity_premium = s.kappa / s.m * int;
thetap = -theta / eta;
sigma_theta = thetap / theta * (psi - eta) * s.sigma;
omega = (1 - eta) / ((1 - eta) + theta * eta);
DS = deposit_spread(psi, omega, eta, s);

% Calculate residuals
% (1) psi = 1 or asset pricing condition
% (2) consumption market-clearing
if psi > 1
    F(1) = log(psi); % force psi to 1
else
    F(1) = log(Q) - log((s.ab - s.ah) / (-sigma_theta * s.sigma + liquidity_premium - DS));
end
F(2) = log(Q) - log(consumption_market_clearing(eta, theta, psi, 'Q', s));
% F(2) = log(Q * (s.r * (1 + (theta - 1) * eta) + s.eps1 / s.eps2)) - ...
%     log((s.ab - s.ah) * psi + s.ah + 1 / s.eps2);

end
