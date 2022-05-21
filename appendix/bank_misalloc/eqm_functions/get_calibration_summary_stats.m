function summary_stats = get_calibration_summary_stats(solution, s)
% This function calculates moments for calibration purposes.
%
% Written by William Chen and Greg Phelan, Jan. 2022

% Ergodic averages
summary_stats.ergodic = get_ergodic_moments(solution, s);

% Moments other than ergodic averages
% summary_stats.other = get_other_moments(solution, s)

end

function out = get_ergodic_moments(sol, s)
% Get ergodic averages of various quantities

    Nd = length(sol.eta_density);

    if Nd == length(sol.eta) - 1 % Then no Infs appeared during ergocalc
        i = 2:(Nd + 1); % See further down in the function

        % Stability
%         indind = sol.psi(2:end) < 1;
%         stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
%         out.stab_frac = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));
%
%         indind = sol.psi(2:end) < .25;
%         stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
%         out.stab_frac_pct25 = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));
%
        indind = sol.psi(2:end) < .5;
        stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
        out.stab_frac_pct50 = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));
%
%         indind = sol.psi(2:end) < .10;
%         stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
%         out.stab_frac_pct10 = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));
%

        % Other calibration quantities
        out.drf = ergodic_mean(sol.drf(i), sol.density, sol.eta_density);
        % out.drd = ergodic_mean(sol.drd(i), sol.density, sol.eta_density);
        out.deposit_spread = ergodic_mean(sol.deposit_spread(i), sol.density, sol.eta_density);
        % out.mu_Ch = ergodic_mean(sol.mu_Ch(i), sol.density, sol.eta_density);
        out.sigma_Ch = ergodic_mean(sol.sigma_Ch(i), sol.density, sol.eta_density);
        % out.mu_K = ergodic_mean(sol.mu_K(i), sol.density, sol.eta_density);
        % out.mu_gdp = ergodic_mean(sol.mu_gdp(i), sol.density, sol.eta_density);
        out.sigma_gdp = ergodic_mean(sol.sigma_gdp(i), sol.density, sol.eta_density);
        out.Sharpe = ergodic_mean(sol.Sharpe(i), sol.density, sol.eta_density);
        % out.risk_premium = ergodic_mean(sol.risk_premium(i), sol.density, sol.eta_density);
        % out.liquidity_premium = ergodic_mean(sol.liquidity_premium(i), sol.density, sol.eta_density);
	    % out.BookLeverageRatio = ergodic_mean(sol.BookLeverageRatio(i), sol.density, sol.eta_density);
	    out.MarketLeverageRatio = ergodic_mean(sol.MarketLeverageRatio(i), sol.density, sol.eta_density);
		out.iota = ergodic_mean(sol.iota(i), sol.density, sol.eta_density);
		% out.mu_iota = ergodic_mean(sol.mu_iota(i), sol.density, sol.eta_density);
		out.sigma_iota = ergodic_mean(sol.sigma_iota(i), sol.density, sol.eta_density);
		out.invst_elasticity = ergodic_mean(sol.invst_elasticity(i), sol.density, sol.eta_density);

%         good_times = sol.psi(i) >= 1;
%         mass_interp = griddedInterpolant(sol.eta_density, sol.density .* good_times);
%         mass = integral(@(x) mass_interp(x), sol.eta_density(1), sol.eta_density(end));
%         good_pi_interp = griddedInterpolant(sol.eta_density, sol.inflation(i) .* good_times .* sol.density ./ mass);
%         out.good_pi = integral(@(x) good_pi_interp(x), sol.eta_density(1), sol.eta_density(end));

        non_zlb = sol.eta_density >= s.etaPUT;
        mass_interp = griddedInterpolant(sol.eta_density, sol.density .* non_zlb);
        mass = integral(@(x) mass_interp(x), sol.eta_density(1), sol.eta_density(end));
        non_zlb_pi_interp = griddedInterpolant(sol.eta_density, sol.inflation(i) .* non_zlb .* sol.density ./ mass);
        out.non_zlb_pi = integral(@(x) non_zlb_pi_interp(x), sol.eta_density(1), sol.eta_density(end));
    else % Handle etas by interpolation
        [~, ia, ib] = intersect(sol.eta, sol.eta_density);
        eta = sol.eta(ia);

        % Stability
%         indind = sol.psi(ia) < 1;
%         stab_frac_interp = griddedInterpolant(eta, indind .* sol.density(ib));
%         out.stab_frac = 100 * (1 - integral(@(x) stab_frac_interp(x), eta(1), eta(end)));

        % Calibration quantities
        out.drf = ergodic_mean(sol.drf(ia), sol.density(ib), eta);
        % out.drd = ergodic_mean(sol.drd(ia), sol.density(ib), eta);
        out.deposit_spread = ergodic_mean(sol.deposit_spread(ia), sol.density(ib), eta);
        % out.mu_Ch = ergodic_mean(sol.mu_Ch(ia), sol.density(ib), eta);
        out.sigma_Ch = ergodic_mean(sol.sigma_Ch(ia), sol.density(ib), eta);
        % out.mu_K = ergodic_mean(sol.mu_K(ia), sol.density(ib), eta);
        out.sigma_gdp = ergodic_mean(sol.sigma_gdp(ia), sol.density(ib), eta);
        % out.mu_gdp = ergodic_mean(sol.mu_gdp(ia), sol.density(ib), eta);
        out.Sharpe = ergodic_mean(sol.Sharpe(ia), sol.density(ib), eta);
        % out.risk_premium = ergodic_mean(sol.risk_premium(ia), sol.density(ib), eta);
        % out.liquidity_premium = ergodic_mean(sol.liquidity_premium(ia), sol.density(ib), eta);
	    % out.BookLeverageRatio = ergodic_mean(sol.BookLeverageRatio(ia), sol.density(ib), eta);
	    out.MarketLeverageRatio = ergodic_mean(sol.MarketLeverageRatio(ia), sol.density(ib), eta);
		out.iota = ergodic_mean(sol.iota(ia), sol.density(ib), eta);
		% out.mu_iota = ergodic_mean(sol.mu_iota(ia), sol.density(ib), eta);
		out.sigma_iota = ergodic_mean(sol.sigma_iota(ia), sol.density(ib), eta);
		out.invst_elasticity = ergodic_mean(sol.invst_elasticity(ia), sol.density(ib), eta);

%         good_times = sol.psi(ia) >= 1;
%         mass_interp = griddedInterpolant(eta, sol.density(ib) .* good_times);
%         mass = integral(@(x) mass_interp(x), eta(1), eta(end));
%         good_pi_interp = griddedInterpolant(eta, sol.inflation(ia) .* good_times .* sol.density(ib) ./ mass);
%         out.good_pi = integral(@(x) good_pi_interp(x), eta(1), eta(end));

        non_zlb = eta >= s.etaPUT;
        mass_interp = griddedInterpolant(eta, sol.density(ib) .* non_zlb);
        mass = integral(@(x) mass_interp(x), eta(1), eta(end));
        non_zlb_pi_interp = griddedInterpolant(eta, sol.inflation(ia) .* non_zlb .* sol.density(ib) ./ mass);
        out.non_zlb_pi = integral(@(x) non_zlb_pi_interp(x), eta(1), eta(end));
    end
end

% function out = get_other_moments(sol, s)
% end
