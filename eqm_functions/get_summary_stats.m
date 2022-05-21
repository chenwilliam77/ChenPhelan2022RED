function summary_stats = get_summary_stats(solution, s)
% This function calculates various summary statistics of interest
% for calibration purposes and presentation of results.
%
% Written by William Chen and Greg Phelan, Jan. 2022

% Uniform averages
% summary_stats.uniform = get_uniform_moments(solution);

if isfield(s, 'max_expV') && s.max_expV
    % Only calculate ex-ante welfare
    Nd = length(solution.eta_density);

    if Nd == length(solution.eta) - 1 % Then no Infs appeared during ergocalc
        i = 2:(Nd + 1); % See further down in the function

        ExpV = ergodic_mean(solution.Welfare(i), solution.density, solution.eta_density);
    else
        [~, ia, ib] = intersect(solution.eta, solution.eta_density);
        eta = solution.eta(ia);
        ExpV = ergodic_mean(solution.Welfare(ia), solution.density(ib), eta);
    end

    ergodic.ExpV = ExpV;
    summary_stats.ergodic = ergodic;
else
    % Ergodic averages
    summary_stats.ergodic = get_ergodic_moments(solution, s);

    % Stochastic steady state
    summary_stats.sss = get_stochastic_steady_state_values(solution, s);
end

end

function out = get_ergodic_moments(sol, s)
% Get ergodic averages of various quantities

    Nd = length(sol.eta_density);

    if Nd == length(sol.eta) - 1 % Then no Infs appeared during ergocalc
        i = 2:(Nd + 1); % See further down in the function

        % Nominal interest rate
        INT = zeros(Nd, 1);
        for n = 1:Nd
            INT(n) = interest(sol.eta_density(n), s);
        end
        interest_interp = griddedInterpolant(sol.eta_density, sol.density .* INT);
        out.interest = integral(@(x) interest_interp(x), sol.eta_density(1), sol.eta_density(end));

        % Stability
        indind = sol.psi(2:end) < 1;
        stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
        out.stab_frac = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));

        indind = sol.psi(2:end) < .25;
        stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
        out.stab_frac_pct25 = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));

        indind = sol.psi(2:end) < .5;
        stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
        out.stab_frac_pct50 = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));

        indind = sol.psi(2:end) < .10;
        stab_frac_interp = griddedInterpolant(sol.eta_density, indind .* sol.density);
        out.stab_frac_pct10 = 100 * (1 - integral(@(x) stab_frac_interp(x), sol.eta_density(1), sol.eta_density(end)));

        % Ex-ante welfare
        out.ExpV      = ergodic_mean(sol.Welfare(i), sol.density, sol.eta_density);
        out.ExpVc     = ergodic_mean(sol.WelfareConsume(i), sol.density, sol.eta_density);
        out.ExpVd     = ergodic_mean(sol.WelfareDeposits(i), sol.density, sol.eta_density);
        out.ExpVpi    = ergodic_mean(sol.WelfarePi(i), sol.density, sol.eta_density);
        out.ExpVc_Vpi = ergodic_mean(sol.WelfareConsume(i) - sol.WelfarePi(i), sol.density, sol.eta_density);

        % Other calibration quantities
        out.drf = ergodic_mean(sol.drf(i), sol.density, sol.eta_density);
        out.drd = ergodic_mean(sol.drd(i), sol.density, sol.eta_density);
        out.deposit_spread = ergodic_mean(sol.deposit_spread(i), sol.density, sol.eta_density);
        out.mu_Ch = ergodic_mean(sol.mu_Ch(i), sol.density, sol.eta_density);
        out.sigma_Ch = ergodic_mean(sol.sigma_Ch(i), sol.density, sol.eta_density);
        out.sigma_Q = ergodic_mean(sol.sigma_Q(i), sol.density, sol.eta_density);
        out.mu_K = ergodic_mean(sol.mu_K(i), sol.density, sol.eta_density);
        out.mu_gdp = ergodic_mean(sol.mu_gdp(i), sol.density, sol.eta_density);
        out.sigma_gdp = ergodic_mean(sol.sigma_gdp(i), sol.density, sol.eta_density);
        out.gdp_capital_value_ratio = ergodic_mean(sol.gdp_capital_value_ratio(i), sol.density, sol.eta_density);
        out.Sharpe = ergodic_mean(sol.Sharpe(i), sol.density, sol.eta_density);
        out.risk_premium = ergodic_mean(sol.risk_premium(i), sol.density, sol.eta_density);
        out.liquidity_premium = ergodic_mean(sol.liquidity_premium(i), sol.density, sol.eta_density);
        out.Leverage = ergodic_mean(sol.Leverage(i), sol.density, sol.eta_density);
	    out.BookLeverageRatio = ergodic_mean(sol.BookLeverageRatio(i), sol.density, sol.eta_density);
	    out.MarketLeverageRatio = ergodic_mean(sol.MarketLeverageRatio(i), sol.density, sol.eta_density);
		out.iota = ergodic_mean(sol.iota(i), sol.density, sol.eta_density);
		out.mu_iota = ergodic_mean(sol.mu_iota(i), sol.density, sol.eta_density);
		out.sigma_iota = ergodic_mean(sol.sigma_iota(i), sol.density, sol.eta_density);
        out.invst_elasticity = ergodic_mean(sol.invst_elasticity(i), sol.density, sol.eta_density);
		out.invst_adjustment_costs = ergodic_mean(sol.invst_adjustment_costs(i), sol.density, sol.eta_density);
		out.invst_capital_value_ratio = ergodic_mean(sol.invst_capital_value_ratio(i), sol.density, sol.eta_density);

        good_times = sol.psi(i) >= 1;
        mass_interp = griddedInterpolant(sol.eta_density, sol.density .* good_times);
        mass = integral(@(x) mass_interp(x), sol.eta_density(1), sol.eta_density(end));
        good_pi_interp = griddedInterpolant(sol.eta_density, sol.inflation(i) .* good_times .* sol.density ./ mass);
        out.good_pi = integral(@(x) good_pi_interp(x), sol.eta_density(1), sol.eta_density(end));
        good_iota_interp = griddedInterpolant(sol.eta_density, sol.iota(i) .* good_times .* sol.density ./ mass);
        out.good_iota = integral(@(x) good_iota_interp(x), sol.eta_density(1), sol.eta_density(end));
        good_invst_cap_ratio_interp = griddedInterpolant(sol.eta_density, sol.invst_capital_value_ratio(i) .* good_times .* sol.density ./ mass);
        out.good_invst_capital_value_ratio = integral(@(x) good_invst_cap_ratio_interp(x), sol.eta_density(1), sol.eta_density(end));
        good_gdp_cap_ratio_interp = griddedInterpolant(sol.eta_density, sol.gdp_capital_value_ratio(i) .* good_times .* sol.density ./ mass);
        out.good_gdp_capital_value_ratio = integral(@(x) good_gdp_cap_ratio_interp(x), sol.eta_density(1), sol.eta_density(end));
        good_booklvg_interp = griddedInterpolant(sol.eta_density, sol.BookLeverageRatio(i) .* good_times .* sol.density ./ mass);
        out.good_BookLeverageRatio = integral(@(x) good_booklvg_interp(x), sol.eta_density(1), sol.eta_density(end));
        good_mktlvg_interp = griddedInterpolant(sol.eta_density, sol.MarketLeverageRatio(i) .* good_times .* sol.density ./ mass);
        out.good_MarketLeverageRatio = integral(@(x) good_mktlvg_interp(x), sol.eta_density(1), sol.eta_density(end));

        non_zlb = sol.eta_density >= s.etaPUT;
        mass_interp = griddedInterpolant(sol.eta_density, sol.density .* non_zlb);
        mass = integral(@(x) mass_interp(x), sol.eta_density(1), sol.eta_density(end));
        non_zlb_pi_interp = griddedInterpolant(sol.eta_density, sol.inflation(i) .* non_zlb .* sol.density ./ mass);
        out.non_zlb_pi = integral(@(x) non_zlb_pi_interp(x), sol.eta_density(1), sol.eta_density(end));
        non_zlb_booklvg_interp = griddedInterpolant(sol.eta_density, sol.BookLeverageRatio(i) .* non_zlb .* sol.density ./ mass);
        out.non_zlb_BookLeverageRatio = integral(@(x) non_zlb_booklvg_interp(x), sol.eta_density(1), sol.eta_density(end));
        non_zlb_mktlvg_interp = griddedInterpolant(sol.eta_density, sol.MarketLeverageRatio(i) .* non_zlb .* sol.density ./ mass);
        out.non_zlb_MarketLeverageRatio = integral(@(x) non_zlb_mktlvg_interp(x), sol.eta_density(1), sol.eta_density(end));
    else % Handle etas by interpolation
        [~, ia, ib] = intersect(sol.eta, sol.eta_density);
        eta = sol.eta(ia);
        Nd = length(eta);

        % Nominal interest rate
        INT = zeros(Nd, 1);
        for n = 1:Nd
            INT(n) = interest(eta(n), s);
        end
        interest_interp = griddedInterpolant(eta, INT .* sol.density(ib));
        out.interest = integral(@(x) interest_interp(x), eta(1), eta(end));

        % Stability
        indind = sol.psi(ia) < 1;
        stab_frac_interp = griddedInterpolant(eta, indind .* sol.density(ib));
        out.stab_frac = 100 * (1 - integral(@(x) stab_frac_interp(x), eta(1), eta(end)));

        indind = sol.psi(ia) < 0.5;
        stab_frac_interp = griddedInterpolant(eta, indind .* sol.density(ib));
        out.stab_frac_pct50 = 100 * (1 - integral(@(x) stab_frac_interp(x), eta(1), eta(end)));

        indind = sol.psi(ia) < 0.25;
        stab_frac_interp = griddedInterpolant(eta, indind .* sol.density(ib));
        out.stab_frac_pct25 = 100 * (1 - integral(@(x) stab_frac_interp(x), eta(1), eta(end)));

        % Ex-ante welfare
        out.ExpV      = ergodic_mean(sol.Welfare(ia), sol.density(ib), eta);
        out.ExpVc     = ergodic_mean(sol.WelfareConsume(ia), sol.density(ib), eta);
        out.ExpVd     = ergodic_mean(sol.WelfareDeposits(ia), sol.density(ib), eta);
        out.ExpVpi    = ergodic_mean(sol.WelfarePi(ia), sol.density(ib), eta);
        out.ExpVc_Vpi = out.ExpVc - out.ExpVpi;

        % Calibration quantities
        out.drf = ergodic_mean(sol.drf(ia), sol.density(ib), eta);
        out.drd = ergodic_mean(sol.drd(ia), sol.density(ib), eta);
        out.deposit_spread = ergodic_mean(sol.deposit_spread(ia), sol.density(ib), eta);
        out.mu_Ch = ergodic_mean(sol.mu_Ch(ia), sol.density(ib), eta);
        out.sigma_Ch = ergodic_mean(sol.sigma_Ch(ia), sol.density(ib), eta);
        out.sigma_Q = ergodic_mean(sol.sigma_Q(ia), sol.density(ib), eta);
        out.mu_K = ergodic_mean(sol.mu_K(ia), sol.density(ib), eta);
        out.sigma_gdp = ergodic_mean(sol.sigma_gdp(ia), sol.density(ib), eta);
        out.mu_gdp = ergodic_mean(sol.mu_gdp(ia), sol.density(ib), eta);
        out.gdp_capital_value_ratio = ergodic_mean(sol.gdp_capital_value_ratio(ia), sol.density(ib), eta);
        out.Sharpe = ergodic_mean(sol.Sharpe(ia), sol.density(ib), eta);
        out.risk_premium = ergodic_mean(sol.risk_premium(ia), sol.density(ib), eta);
        out.liquidity_premium = ergodic_mean(sol.liquidity_premium(ia), sol.density(ib), eta);
        out.Leverage = ergodic_mean(sol.Leverage(ia), sol.density(ib), eta);
	    out.BookLeverageRatio = ergodic_mean(sol.BookLeverageRatio(ia), sol.density(ib), eta);
	    out.MarketLeverageRatio = ergodic_mean(sol.MarketLeverageRatio(ia), sol.density(ib), eta);
		out.iota = ergodic_mean(sol.iota(ia), sol.density(ib), eta);
		out.mu_iota = ergodic_mean(sol.mu_iota(ia), sol.density(ib), eta);
		out.sigma_iota = ergodic_mean(sol.sigma_iota(ia), sol.density(ib), eta);
        out.invst_elasticity = ergodic_mean(sol.invst_elasticity(ia), sol.density(ib), eta);
		out.invst_adjustment_costs = ergodic_mean(sol.invst_adjustment_costs(ia), sol.density(ib), eta);
		out.invst_capital_value_ratio = ergodic_mean(sol.invst_capital_value_ratio(ia), sol.density(ib), eta);

        good_times = sol.psi(ia) >= 1;
        mass_interp = griddedInterpolant(eta, sol.density(ib) .* good_times);
        mass = integral(@(x) mass_interp(x), eta(1), eta(end));
        good_pi_interp = griddedInterpolant(eta, sol.inflation(ia) .* good_times .* sol.density(ib) ./ mass);
        out.good_pi = integral(@(x) good_pi_interp(x), eta(1), eta(end));
        good_iota_interp = griddedInterpolant(eta, sol.iota(ia) .* good_times .* sol.density(ib) ./ mass);
        out.good_iota = integral(@(x) good_iota_interp(x), eta(1), eta(end));
        good_invst_cap_ratio_interp = griddedInterpolant(eta, sol.invst_capital_value_ratio(ia) .* good_times .* sol.density(ib) ./ mass);
        out.good_invst_capital_value_ratio = integral(@(x) good_invst_cap_ratio_interp(x), eta(1), eta(end));
        good_gdp_cap_ratio_interp = griddedInterpolant(eta, sol.gdp_capital_value_ratio(ia) .* good_times .* sol.density(ib) ./ mass);
        out.good_gdp_capital_value_ratio = integral(@(x) good_gdp_cap_ratio_interp(x), eta(1), eta(end));
        good_booklvg_interp = griddedInterpolant(eta, sol.BookLeverageRatio(ia) .* good_times .* sol.density(ib) ./ mass);
        out.good_BookLeverageRatio = integral(@(x) good_booklvg_interp(x), eta(1), eta(end));
        good_mktlvg_interp = griddedInterpolant(eta, sol.MarketLeverageRatio(ia) .* good_times .* sol.density(ib) ./ mass);
        out.good_MarketLeverageRatio = integral(@(x) good_mktlvg_interp(x), eta(1), eta(end));

        non_zlb = eta >= s.etaPUT;
        mass_interp = griddedInterpolant(eta, sol.density(ib) .* non_zlb);
        mass = integral(@(x) mass_interp(x), eta(1), eta(end));
        non_zlb_pi_interp = griddedInterpolant(eta, sol.inflation(ia) .* non_zlb .* sol.density(ib) ./ mass);
        out.non_zlb_pi = integral(@(x) non_zlb_pi_interp(x), eta(1), eta(end));
        non_zlb_booklvg_interp = griddedInterpolant(eta, sol.BookLeverageRatio(ia) .* non_zlb .* sol.density(ib) ./ mass);
        out.non_zlb_BookLeverageRatio = integral(@(x) non_zlb_booklvg_interp(x), eta(1), eta(end));
        non_zlb_mktlvg_interp = griddedInterpolant(eta, sol.MarketLeverageRatio(ia) .* non_zlb .* sol.density(ib) ./ mass);
        out.non_zlb_MarketLeverageRatio = integral(@(x) non_zlb_mktlvg_interp(x), eta(1), eta(end));
    end
end

function out = get_stochastic_steady_state_values(sol, s)

    sss_il = false(length(sol.eta), 1); % Left index
    sss_ir = false(length(sol.eta), 1); % right index
    w_il   = zeros(length(sol.eta), 1);  % weight on left index
    w_ir   = zeros(length(sol.eta), 1);  % weight on right index
    if all(sol.mu_eta > 0) % then right boundary is SSS
        w_ir(end) = 1;
        sss_ir(end) = true;
    else
        % interior SSS (it is assumed that mu_eta is not always negative,
        %               in which case parameters are poorly chosen)
	    for i = 1:(length(sol.eta) - 1)
	        % Are there opposite signs?
	        if sol.mu_eta(i) * sol.mu_eta(i + 1) < 0
	            sss_il(i) = true;
	            sss_ir(i + 1) = true;
	            w_il(i) = abs(sol.mu_eta(i)) / (abs(sol.mu_eta(i)) + abs(sol.mu_eta(i + 1)));
	            w_ir(i + 1) = 1 - w_il(i);
	        end
	    end
    end

    out.interest = w_il(sss_il) .* sol.interestvec(sss_il) + w_ir(sss_ir) .* sol.interestvec(sss_ir);
    out.drf = w_il(sss_il) .* sol.drf(sss_il) + w_ir(sss_ir) .* sol.drf(sss_ir);
    out.mu_Ch = w_il(sss_il) .* sol.mu_Ch(sss_il) + w_ir(sss_ir) .* sol.mu_Ch(sss_ir);
    out.sigma_Ch = w_il(sss_il) .* sol.sigma_Ch(sss_il) + w_ir(sss_ir) .* sol.sigma_Ch(sss_ir);
    out.sigma_Q = w_il(sss_il) .* sol.sigma_Q(sss_il) + w_ir(sss_ir) .* sol.sigma_Q(sss_ir);
    out.mu_K = w_il(sss_il) .* sol.mu_K(sss_il) + w_ir(sss_ir) .* sol.mu_K(sss_ir);
    out.gdp = w_il(sss_il) .* sol.gdp(sss_il) + w_ir(sss_ir) .* sol.gdp(sss_ir);
    out.mu_gdp = w_il(sss_il) .* sol.mu_gdp(sss_il) + w_ir(sss_ir) .* sol.mu_gdp(sss_ir);
    out.sigma_gdp = w_il(sss_il) .* sol.sigma_gdp(sss_il) + w_ir(sss_ir) .* sol.sigma_gdp(sss_ir);
    out.Sharpe = w_il(sss_il) .* sol.Sharpe(sss_il) + w_ir(sss_ir) .* sol.Sharpe(sss_ir);
    out.Leverage = w_il(sss_il) .* sol.Leverage(sss_il) + w_ir(sss_ir) .* sol.Leverage(sss_ir);
    out.BookLeverageRatio = w_il(sss_il) .* sol.BookLeverageRatio(sss_il) + w_ir(sss_ir) .* sol.BookLeverageRatio(sss_ir);
    out.MarketLeverageRatio = w_il(sss_il) .* sol.MarketLeverageRatio(sss_il) + w_ir(sss_ir) .* sol.MarketLeverageRatio(sss_ir);
    out.Welfare = w_il(sss_il) .* sol.Welfare(sss_il) + w_ir(sss_ir) .* sol.Welfare(sss_ir);
    out.WelfareConsume = w_il(sss_il) .* sol.WelfareConsume(sss_il) + w_ir(sss_ir) .* sol.WelfareConsume(sss_ir);
    out.WelfareDeposits = w_il(sss_il) .* sol.WelfareDeposits(sss_il) + w_ir(sss_ir) .* sol.WelfareDeposits(sss_ir);
    out.WelfarePi = w_il(sss_il) .* sol.WelfarePi(sss_il) + w_ir(sss_ir) .* sol.WelfarePi(sss_ir);
    out.eta = w_il(sss_il) .* sol.eta(sss_il) + w_ir(sss_ir) .* sol.eta(sss_ir);
    out.cdf = interp1(sol.eta_density, sol.cdf, out.eta, s.interp);
end

% function out = uniform_mean(x, y)
%     out = dot(x, [y(1); diff(y)]) ./ y(end);
% end
%
% function out = get_uniform_moments(sol)
%
%     % Nominal interest rate
%     out.interest = uniform_mean(sol.interestvec, sol.eta);
%
%     % Calibration quantities
%     out.drf = uniform_mean(sol.drf, sol.eta);
%     out.mu_Ch = uniform_mean(sol.mu_Ch, sol.eta);
%     out.sigma_Ch = uniform_mean(sol.sigma_Ch, sol.eta);
%     out.sigma_Q = uniform_mean(sol.sigma_Q, sol.eta);
%     out.mu_K = uniform_mean(sol.mu_K, sol.eta);
%     out.mu_gdp = uniform_mean(sol.mu_gdp, sol.eta);
%     out.sigma_gdp = uniform_mean(sol.sigma_gdp, sol.eta);
%     out.Sharpe = uniform_mean(sol.Sharpe, sol.eta);
%     out.good_pi = uniform_mean(sol.inflation, sol.eta);
%     out.Leverage = uniform_mean(sol.Leverage, sol.eta);
%     out.BookLeverageRatio = uniform_mean(sol.BookLeverageRatio, sol.eta);
%     out.MarketLeverageRatio = uniform_mean(sol.MarketLeverageRatio, sol.eta);
% 	out.iota = uniform_mean(sol.iota, sol.eta);
% end
