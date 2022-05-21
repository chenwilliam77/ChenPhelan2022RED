  function plot_results(solution, s)
    % Set up
    color       = s.color;
    N           = length(solution.eta);
    etastar     = solution.eta(end);
    start       = s.start;

    % Plot
    figure(1); hold on
        ylim([0 50]);
        plot(solution.eta, solution.Leverage, 'color', color, 'linestyle', s.line, 'linewidth', s.width);
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b', 'Interpreter', 'latex');
        ylabel('Leverage', 'FontSize', s.fontsize, 'FontWeight', 'b');
        if isfield(s, 'use_title')
            if s.use_title
		        title('Aggregate Bank Leverage','FontSize', s.fontsize);
			end
		end
        xlim([s.start 1.05 * etastar]);

    figure(2); hold on
        plot(solution.eta, solution.mu_eta .* solution.eta, 'color',color,'linestyle',s.line,'linewidth',s.width);
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b', 'Interpreter', 'latex');
        ylabel('$\eta \mu_\eta$ ','fontsize',s.fontsize,'FontWeight','b');
        xlim([start 1.05 * etastar]);
        if isfield(s, 'use_title')
            if s.use_title
		        title('Drift of bank equity','FontSize',s.fontsize);
            end
		end

    figure(22); hold on
        plot(solution.eta, solution.sigma_eta .* solution.eta, 'color',color,'linestyle',s.line,'linewidth',s.width);
        ylabel('$\eta \sigma_\eta$ ','fontsize',s.fontsize,'FontWeight','b');
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b', 'Interpreter', 'latex');
        xlim([start 1.05 * etastar]) ;
        % ylim([0 1.05*max(solution.sigma_eta .* solution.eta)]);
        if isfield(s, 'use_title')
            if s.use_title
                title('Volatility of bank equity','FontSize',s.fontsize);
            end
        end

    figure(23); hold on
        plot(solution.eta, solution.psi, 'color',color,'linestyle',s.line,'linewidth',s.width);
        ylabel('$\psi$ ','fontsize',s.fontsize,'FontWeight','b');
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b', 'Interpreter', 'latex');
        xlim([start 1.05 * etastar]);
        if isfield(s, 'use_title')
            if s.use_title
                title('Share of Capital Held by Banks','FontSize',s.fontsize);
            end
        end

    figure(5); hold on
        plot(solution.eta, solution.flow_util, 'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('Flow Utility','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Flow Utility','FontSize',s.fontsize);
            end
        end

        xlim([start 1.05*etastar])

    figure(6); hold on
        plot(solution.eta, solution.omega, 'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\omega$','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Share of Household Wealth in Non-Bank Net Worth', 'FontSize', s.fontsize, 'FontWeight', 'b')
            end
        end
        xlim([start 1.05*etastar])

    figure(7); hold on
        plot(solution.eta, solution.sigma_Ch, 'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\sigma_{Ch}$','FontSize',s.fontsize,'FontWeight','b')
        title('Volatility of Household Consumption', 'FontSize', s.fontsize, 'FontWeight', 'b');
        xlim([start 1.05*etastar])

    figure(8); hold on
        plot(solution.eta, solution.sigma_theta, 'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\sigma_{\theta}$','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Volatility of Marginal Value of Bank Equity', 'FontSize', s.fontsize, 'FontWeight', 'b');
            end
        end
        xlim([start 1.05*etastar])

    figure(9); hold on
	    plot(solution.eta, solution.Q .* solution.sigma_Q, 'color',color,'linestyle',s.line,'linewidth',s.width);
	    ylabel('$Q \sigma_Q $','fontsize',s.fontsize,'FontWeight','b');
	    xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b', 'Interpreter', 'latex');
	    xlim([s.start 1.05 * etastar]);
	    % ylim([0 1.05*max(solution.sigma_Q .* solution.Q)]);
	    if isfield(s, 'use_title')
	        if s.use_title
	            title('Volatility of price Q','FontSize',s.fontsize);
	        end
	    end

    figure(10); hold on
        plot(solution.eta, solution.Q, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Price $Q$','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$Q$','FontSize',s.fontsize,'FontWeight','b')

    figure(11); hold on
        plot(solution.eta, solution.drf, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Risk-free Interest Rate','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$dr_f$','FontSize',s.fontsize,'FontWeight','b')

    figure(12); hold on
        plot(solution.eta, solution.Welfare, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Household Welfare','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$V(\eta)$','FontSize',s.fontsize,'FontWeight','b')

    figure(14); hold on
        plot(solution.eta, solution.theta, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Marginal Value of Bank Equity','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\theta(\eta)$','FontSize',s.fontsize,'FontWeight','b')

    figure(16)
        plot(solution.eta, solution.mu_K,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\mu_{K}$','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Capital Growth Rate');
            end
        end
        xlim([start 1.05*etastar])

    figure(18)
        plot(solution.eta, solution.iota, 'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\iota$','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Investment Rate');
            end
        end
        xlim([start 1.05*etastar])

    figure(19)
        plot(solution.eta, solution.mu_theta,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$\mu_{\theta}$','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('$\theta$ Growth Rate');
            end
        end
        xlim([start 1.05*etastar])

    figure(20); hold on
        plot(solution.eta, solution.BookLeverageRatio, 'color', color, 'linestyle', s.line, 'linewidth', s.width);
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b', 'Interpreter', 'latex');
        ylabel('Book Leverage', 'FontSize', s.fontsize, 'FontWeight', 'b');
        if isfield(s, 'use_title')
            if s.use_title
		        title('Aggregate Bank Book Leverage Ratio','FontSize', s.fontsize);
			end
		end
        xlim([s.start 1.05 * etastar]);

    figure(21); hold on
        plot(solution.eta, solution.MarketLeverageRatio, 'color', color, 'linestyle', s.line, 'linewidth', s.width);
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share', 'FontSize', s.fontsize, 'FontWeight', 'b', 'Interpreter', 'latex');
        ylabel('Market Leverage', 'FontSize', s.fontsize, 'FontWeight', 'b');
        if isfield(s, 'use_title')
            if s.use_title
		        title('Aggregate Bank Market Leverage Ratio','FontSize', s.fontsize);
			end
		end
        xlim([s.start 1.05 * etastar]);

    figure(24)
        plot(solution.eta, solution.Ch, 'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta$, level of bank equity','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$C_h$','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Household Consumption');
            end
        end
        xlim([start 1.05*etastar])

    figure(25)
        plot(solution.eta,solution.Sharpe,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('SR','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Sharpe Ratio')
            end
        end
        xlim([start 1.05 * etastar])

    figure(26)
        plot(solution.eta,solution.risk_premium + solution.liquidity_premium,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('Excess Returns','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title("Banks' Excess Returns")
            end
        end
        xlim([start 1.05 * etastar])


    figure(32),
        if s.nomruletype == 1
            ind_put = find(solution.eta >= s.etaPUT, 1);
            ind_law = find(solution.eta >= s.etaLAW, 1);
            ploteta = [solution.eta(1:ind_put - 1); NaN; solution.eta(ind_put:ind_law - 1); NaN; solution.eta(ind_law:end)];
            plotint = [solution.interestvec(1:ind_put - 1); NaN; solution.interestvec(ind_put:ind_law - 1); NaN; solution.interestvec(ind_law:end)];
            plot(ploteta,100 * plotint,'color',color,'linestyle', s.line, 'linewidth',s.width); hold on
        else
            plot(solution.eta,100 * solution.interestvec,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        end
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('\%','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Nominal Interest Rate')
            end
        end
        xlim([start 1.05*etastar])

    figure(33),
        plot(solution.eta,100 * solution.inflation,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('annualized \%','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Inflation Rate');
            end
        end
        xlim([start 1.05*etastar])

    figure(35)
        plot(solution.eta,solution.risk_premium,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        plot(solution.eta,solution.liquidity_premium,'color',color,'linestyle',':','linewidth',s.width); hold on

        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('risk/liquidity premium','FontSize',s.fontsize,'FontWeight','b')
        if isfield(s, 'use_title')
            if s.use_title
                title('Risk and Liquidity Premia')
            end
        end
        legend('SP','LP','location','best')

    figure(40)
        plot(solution.eta_density, solution.density,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        xlim([s.start 1.05*etastar])
        ylabel('$f(\eta)$','fontsize',s.fontsize,'FontWeight','b');
        if isfield(s, 'use_title')
            if s.use_title
                title('Stationary Density of Bank Equity Levels','FontSize',s.fontsize)
            end
        end
        ylim([0 max(solution.density) * 1.05])

        figure(41)
        plot(solution.eta_density ./ solution.eta_density(end), solution.eta_density(end) * solution.density,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta/\eta^*$, Normalized','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$f(\eta / \eta^*)$','fontsize',s.fontsize,'FontWeight','b');
        if isfield(s, 'use_title')
            if s.use_title
                title('Stationary Density of Bank Equity Levels','FontSize',s.fontsize)
            end
        end

    figure(42)
        plot(solution.eta_density, solution.cdf,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        xlim([s.start 1.05*etastar])
        ylabel('$F(\eta)$','fontsize',s.fontsize,'FontWeight','b');
        if isfield(s, 'use_title')
            if s.use_title
                title('Cumulative Density','FontSize',s.fontsize)
            end
        end
        ylim([0 1])

    figure(43)
        plot(solution.eta_density./ solution.eta_density(end), solution.cdf,'color',color,'linestyle',s.line,'linewidth',s.width); hold on
        xlabel('$\eta/\eta^*$, Normalized','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$F(\eta / \eta^*)$','fontsize',s.fontsize,'FontWeight','b');
        if isfield(s, 'use_title')
            if s.use_title
                title('Cumulative Density','FontSize',s.fontsize)
            end
        end
        ylim([0 1])

    figure(51); hold on
        plot(solution.eta, solution.drd, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Deposit Interest Rate','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$dr_d$','FontSize',s.fontsize,'FontWeight','b')

    figure(52); hold on
        plot(solution.eta, solution.deposit_spread, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Deposit Spread','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('Deposit Spread','FontSize',s.fontsize,'FontWeight','b')

    figure(54); hold on
        plot(solution.eta, solution.convenience_yield, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Convenience Yield','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$v(x_{hd})$','FontSize',s.fontsize,'FontWeight','b')

    figure(55); hold on
        plot(solution.eta, solution.WelfareConsume - solution.WelfarePi, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Welfare Excl. Convenience Yield','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$V_c - V_\pi$','FontSize',s.fontsize,'FontWeight','b')

    figure(56); hold on
        plot(solution.eta, solution.WelfareConsume, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Welfare from Consumption','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$V_c$','FontSize',s.fontsize,'FontWeight','b')

    figure(57); hold on
        plot(solution.eta, solution.WelfareDeposits, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Welfare from Convenience Yield','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$V_d$','FontSize',s.fontsize,'FontWeight','b')

    figure(58); hold on
        plot(solution.eta, solution.WelfarePi, 'color',color,'linestyle',s.line,'linewidth',s.width);
        if isfield(s, 'use_title')
            if s.use_title
                title('Welfare from Inflation Costs','FontSize',s.fontsize)
            end
        end
        xlim([start 1.05 * etastar])
        xlabel('$\eta = \frac{N_b}{QK}$, bank wealth share','FontSize',s.fontsize,'FontWeight','b')
        ylabel('$V_\pi$','FontSize',s.fontsize,'FontWeight','b')

end
