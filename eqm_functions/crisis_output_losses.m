function avg_gdp_loss = crisis_output_losses(sol, T, N, s)
% avg_gdp_loss = crisis_output_losses(sol, T, N, s)
%
% Calculate the average output loss in GDP
% over T years after a crisis starts
% using the Euler-Maruyama algorithm.
% The output loss is relative to staying in the stochastic steady state.
%
% sol: solution struct
%
% T: number of years
%
% N: number of simulation iterations
%
% s: settings struct
%
% Written by William Chen and Gregory Phelan, Jan. 2022

% Calculate total number of periods
dt = 1 / 12; % dt = 1 month
nT = T / dt; % convert years to months

% Find crisis initial eta
crisis_i = find(sol.psi >= 0.5, 1);
crisis_eta = sol.eta(crisis_i);

% Initialize storage containers
dW = normrnd(0, sqrt(dt), nT, N); % Draw Brownian motion shocks with discretization dt
GDP_stock = NaN(N, 1); % container for GDP stock over T years
GDP_sss_stock = NaN(N, 1); % container for GDP stock over T years in the stochastic steady state

% Initialize interpolation functions
interp_eta_mul_mu_eta = griddedInterpolant(sol.eta, sol.mu_eta .* sol.eta, s.interp);
interp_eta_mul_sigma_eta = griddedInterpolant(sol.eta, sol.sigma_eta .* sol.eta, s.interp);
interp_gdp = griddedInterpolant(sol.eta, sol.gdp, s.interp);
interp_mu_K = griddedInterpolant(sol.eta, sol.mu_K, s.interp);
eta_sss = find(sol.mu_eta < 0, 1);
mu_K_sss = interp_mu_K((sol.eta(eta_sss) + sol.eta(eta_sss - 1)) / 2);

% Perform Euler-Maruyama recursion
etastar = sol.eta(end); % needed to implement reflecting boundary at eta*

GDP_flow = NaN(nT, 1); % container for GDP flow over T years w/in each simulation iteration
K_sss_stock = NaN(nT, 1); % container for stock of capital over time when fixed in the stochastic steady state
for sim_iter = 1:N
    eta = crisis_eta; % initialize eta location
    K = 1; % initialize K at 1
    K_sss = 1; % for stochastic steady state results
    for t = 1:nT
        eta_mul_mu_eta = interp_eta_mul_mu_eta(eta);
        eta_mul_sigma_eta = interp_eta_mul_sigma_eta(eta);
        mu_K_K = interp_mu_K(eta) * K;
        eta = euler_maruyama(eta, eta_mul_mu_eta, eta_mul_sigma_eta, dt, dW(t, sim_iter));
        K = euler_maruyama(K, mu_K_K, s.sigma * K, dt, dW(t, sim_iter));
        K_sss = euler_maruyama(K_sss, mu_K_sss, s.sigma * K_sss, dt, dW(t, sim_iter));

        % Handle boundary violations by assuming
        % distance traveled across boundary is how much it is reflected
        if eta < s.start
            eta = s.start - (eta - s.start);
        elseif eta > etastar
            eta = etastar - (eta - etastar);
        end

        % calculate GDP flow
        GDP_flow(t) = interp_gdp(eta) * K;
        K_sss_stock(t) = K_sss;
    end

    % Compute GDP stock over T years
    GDP_stock(sim_iter) = dt * sum(GDP_flow);
    GDP_sss_stock(sim_iter) = dt * sum(s.ab * K_sss_stock);
end

% Calculate GDP losses relative to staying in normal times
avg_gdp_loss = mean(GDP_stock ./ GDP_sss_stock) - 1;
end

function y_next = euler_maruyama(y_prev, mu, sigma, dt, dW)
% Implements one step of the Euler iteration scheme

y_next = y_prev + mu * dt + sigma * dW;
end

% speed ups from parallelizing this step is likely minimal, so the code is
% below but commented out
%
% if isfield(s, 'parallel') && s.parallel
%     parfor sim_iter = 1:N
%         eta = crisis_eta; % initialize eta location
%         GDP_flow = NaN(nT, 1); % container for GDP flow over T years w/in each simulation iteration, declare here b/c parallelization
%         for t = 1:nT
%             eta_mul_mu_eta = interp_eta_mul_mu_eta(eta);
%             eta_mul_sigma_eta = interp_eta_mul_sigma_eta(eta);
%             eta = euler_maruyama(eta, eta_mul_mu_eta, eta_mul_sigma_eta, dt, dW(t, sim_iter));
%
%             % Handle boundary violations by assuming
%             % distance traveled across boundary is how much it is reflected
%             if eta < s.start
%                 eta = s.start - (eta - s.start);
%             elseif eta > etastar
%                 eta = etastar - (eta - etastar);
%             end
%
%             % calculate GDP flow
%             GDP_flow(t) = interp_gdp(eta);
%         end
%
%         % Compute GDP stock over T years
%         GDP_stock(sim_iter) = dt * sum(GDP_flow);
%     end
% else
