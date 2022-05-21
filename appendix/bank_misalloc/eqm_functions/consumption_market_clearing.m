function out = consumption_market_clearing(eta, theta, x, to_calc, s)
% Implements the consumption market-clearing condition
% and returns either Q or psi.
% Version with moral hazard extension
%
% eta:     state variable
%
% theta:   marginal value of bank equity
%
% x:       either Q or psi
%
% to_calc: string specifying what to calculate.
%          Must be one of the following strings -
%          ['Q', 'psi']
%
% s:       settings struct
%
% Written by William Chen and Greg Phelan, Jan. 2022

if isfield(s, 'Phi_form')
    Phi_form = s.Phi_form;
else
    Phi_form = 1;
end

if Phi_form == 1
    ab = productivity_fnct(eta, 0, s);
    deprec_adjust = 1 - s.eps2 * s.delta; % investing iota => no adjustment costs
    denom = s.r * (1 - eta + theta .* eta) + s.eps1 / s.eps2;
    if strcmp(to_calc, 'Q')
        out = (s.ah + (ab - s.ah) * x + deprec_adjust / s.eps2) ./ denom;
    elseif strcmp(to_calc, 'psi')
        out = (x .* denom - s.ah - deprec_adjust / s.eps2) / (ab - s.ah);
    else
        error(['Input to_calc must be either Q or psi; it cannot be ' to_calc]);
    end
elseif Phi_form == 2
    ab = productivity_fnct(eta, 0, s);
    deprec_adjust = 1; % investing iota => no adjustment costs
    denom = s.r * (1 + (theta - 1) .* eta) + s.eps1 / s.eps2;
    if strcmp(to_calc, 'Q')
        out = (s.ah + (ab - s.ah) * x + deprec_adjust / s.eps2) ./ denom;
    elseif strcmp(to_calc, 'psi')
        out = (x .* denom - s.ah - deprec_adjust / s.eps2) / (ab - s.ah);
    else
        error(['Input to_calc must be either Q or psi; it cannot be ' to_calc]);
    end
end

end
