function [opt_i, exitflag] = compute_first_best_rate(s)
% out = compute_first_best_rate(s)
%
% Calculate the first best interest rate
%
% Written by William Chen and Greg Phelan, Feb. 2022

% Calculate optimal i if allocating all capital to banks 
gamma = s.kappa / s.m;
Q1 = consumption_market_clearing(0, 1, 1, 'Q', s);
Phi1 = investment_fnct(Q1, 'Phi', s);
i1 = (s.r + Phi1 - s.delta - s.sigma^2 + s.pi_target) / (1 + gamma);

% Check if it is implementable to allocate all capital to banks
if i1 <= (s.ab - s.ah) / (Q1 * gamma)
    % Implementable!
    opt_i = i1;
    exitflag = 1;
else
    % Numerically search for maximum
    opt = optimoptions('fmincon', 'Display', 'none');
    psi_init = .5;
    [opt_psi, ~, exitflag, ~] = fmincon(@(psi) obj_fnct(psi, s), psi_init, ...
        [], [], [], [], -.1, 1, [], opt);

    % Check if the optimal psi violates short sale constraint
    if opt_psi <= 0
        % Optimal to allocate all capital to HH
        Q = consumption_market_clearing(0, 1, 0, 'Q', s);
        Phi = investment_fnct(Q, 'Phi', s);
        opt_i = s.r + Phi - s.delta - s.sigma^2 + s.pi_target;
    else
        Q = consumption_market_clearing(0, 1, opt_psi, 'Q', s);
        opt_i = (s.ab - s.ah) / (gamma * Q);
    end
end

% Enforce ZLB
if opt_i < 0
    opt_i = 0;
end

end

function F = obj_fnct(psi, s)
gamma = s.kappa / s.m;

% Calculate Q and i given psi
Q = consumption_market_clearing(0, 1, psi, 'Q', s);
i = (s.ab - s.ah) / (gamma * Q);

% Calculate objective function
Phi = investment_fnct(Q, 'Phi', s);
rf = s.r + Phi - s.delta - s.sigma^2 - psi * gamma * i;

F = log(s.r * Q) - s.lambda / 2 * (i - rf - s.pi_target)^2 + Phi / s.r;
F = -F; % take negative b/c minimizer algorithm used

end