function out = ergodic_mean(x, density, eta)
% Compute the ergodic mean (integrate w.r.t. stationary density)
%
% Written by William Chen and Greg Phelan, Jan. 2022

    % mass = density' * diff(eta); % Naive Riemann sum approach
    % out  = (x .* density)' * diff(eta) / mass;
    mass_interp = griddedInterpolant(eta, density); % integrating using interpolation and GK quadrature
    mass = integral(@(x) mass_interp(x), eta(1), eta(end));
    mean_interp = griddedInterpolant(eta, density .* x ./ mass);
    out  = integral(@(x) mean_interp(x), eta(1), eta(end));
end
