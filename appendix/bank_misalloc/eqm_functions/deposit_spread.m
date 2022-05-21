function out = deposit_spread(psi, omega, eta, s)
% Computes the deposit spread due to
% the convenience yield on bank deposits
% in households' utility function.
%
% Written by William Chen and Greg Phelan, Jan. 2022

x_hd = (psi - eta) ./ (1 - eta);
out = s.r ./ omega .* convenience_yield(x_hd, 1, s);

end
