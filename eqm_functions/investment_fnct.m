function out = investment_fnct(Q, to_calc, s)
% Implements the internal investment function
% for capital growth. This function computes more than
% Phi(iota), however, and is intended to centralize
% the calculation of the investment-related quantities.
%
% The function implemented here is
%
% Phi(iota) = (eps1 / eps2) * log((eps2 * iota + exp(eps2 * delta) - eps2 * delta) / eps1)
%
% To recover the form used in the paper, set eps1 = 1.
%
% Q:       capital price
%
% to_calc: string specifying what to calculate.
%          Must be one of the following strings -
%          ['Phi', 'iota', 'adjustment_costs', 'invst_elasticity']
%
% s:       settings struct
%
% Written by William Chen and Greg Phelan, Jan. 2022

% Check if exogenous TFP growth is a parameter
if isfield(s, 'mu_a')
   mu_a = s.mu_a;
else % Set to zero otherwise
   mu_a = 0;
end

if isfield(s, 'Phi_form')
    Phi_form = s.Phi_form;
else
    Phi_form = 1;
end

if Phi_form == 1 % Phi(iota) = delta + log(eps * iota + 1 - eps * delta) / eps
    if strcmp(to_calc, 'Phi')
        out = mu_a + s.delta + (s.eps1 / s.eps2) * log(Q);
    elseif strcmp(to_calc, 'iota')
        deprec_adjust = 1 - s.eps2 * s.delta;
        out = (s.eps1 * Q - deprec_adjust) / s.eps2;
    elseif strcmp(to_calc, 'iotap')
        out = s.eps1 / s.eps2;
    elseif strcmp(to_calc, 'adjustment_costs')
        deprec_adjust = 1 - s.eps2 * s.delta;
        Phi = mu_a + s.delta + (s.eps1 / s.eps2) * log(Q);
        out = (1 / s.eps2) * (exp(s.eps2 * (Phi - s.delta)) - deprec_adjust); % invest this amount to get Phi
    elseif strcmp(to_calc, 'invst_elasticity')
        deprec_adjust = 1 - s.eps2 * s.delta;
        iota = (s.eps1 * Q - deprec_adjust) / s.eps2;
        out = s.eps2 * iota ./ (s.eps2 * iota + deprec_adjust);
    end
elseif Phi_form == 2
    if strcmp(to_calc, 'Phi')
        out = mu_a + (s.eps1 / s.eps2) * log(Q);
    elseif strcmp(to_calc, 'iota')
        deprec_adjust = 1;
        out = (s.eps1 * Q - deprec_adjust) / s.eps2;
    elseif strcmp(to_calc, 'iotap')
        out = s.eps1 / s.eps2;
    elseif strcmp(to_calc, 'adjustment_costs')
        deprec_adjust = 1;
        Phi = mu_a + (s.eps1 / s.eps2) * log(Q);
        out = (1 / s.eps2) * (exp(s.eps2 * Phi) - deprec_adjust); % invest this amount to get Phi
    elseif strcmp(to_calc, 'invst_elasticity')
        deprec_adjust = 1;
        iota = (s.eps1 * Q - deprec_adjust) / s.eps2;
        out = s.eps2 * iota ./ (s.eps2 * iota + deprec_adjust);
    end
end

end
