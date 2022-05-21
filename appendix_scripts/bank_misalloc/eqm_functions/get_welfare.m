function [V, H, Vc, Vd, Vpi] = get_welfare(orig_eta, fout, Dyn, other_dyn, s)
% Calculate welfare using a finite difference scheme.
% Outputs are
% 1. V   (total welfare)
% 2. H   (defined in text and below)
% 3. Vc  (only flow consumption utility)
% 4. Vd  (only flow convenience yield on deposits)
% 5. Vpi (only flow inflation costs)
%
% Welfare is computed as follows:
%
% V(K, eta) = log(K) / r + H(eta) (where H = \log(Q(1 + (\theta - 1)\eta)) / r + h)
%
% and we subsume the dependence on K by normalizing K = 1 so that
%
% V(K, eta) = V(eta) = H(eta)
%
% Under alternative definitions of H, e.g. H = log(Q) / r + h,
% V(eta) may not necessarily coincide with H(eta).
% We found that our definition of H was the most numerically
% accurate. For example, we verified that in economies where
% V(eta) should be constant, defining H as above
% always resulted in an approximately constant V(eta).
%
% Written by William Chen and Greg Phelan, Jan. 2022

% Grab values
orig_sigma_eta       = Dyn(:, 2);
orig_mu_eta          = Dyn(:, 4);
orig_mu_theta        = Dyn(:, 7);
orig_sigma_theta     = Dyn(:, 8);
orig_theta           = fout(:, 1);
orig_Q               = fout(:, 3);
orig_eta_mu_eta      = orig_eta .* orig_mu_eta;
orig_eta_sigma_etasq = (orig_eta .* orig_sigma_eta).^2;
orig_x_hd            = (Dyn(:, 1) - orig_eta) ./ (1 - orig_eta);

% Interpolate to uniform grid so the centered difference
% approximation is properly constructed
eta             = linspace(orig_eta(1), orig_eta(end), s.N_welfare)';
mu_eta          = interp1(orig_eta, orig_mu_eta,          eta, s.interp);
sigma_eta       = interp1(orig_eta, orig_sigma_eta,       eta, s.interp);
mu_theta        = interp1(orig_eta, orig_mu_theta,        eta, s.interp);
sigma_theta     = interp1(orig_eta, orig_sigma_theta,     eta, s.interp);
theta           = interp1(orig_eta, orig_theta,           eta, s.interp);
Q               = interp1(orig_eta, orig_Q,               eta, s.interp);
eta_mu_eta      = interp1(orig_eta, orig_eta_mu_eta,      eta, s.interp);
eta_sigma_etasq = interp1(orig_eta, orig_eta_sigma_etasq, eta, s.interp);
x_hd            = interp1(orig_eta, orig_x_hd,            eta, s.interp);
inflation       = interp1(orig_eta, other_dyn.inflation,  eta, s.interp);

% Set up finite difference operators via upwind scheme
% with reflecting boundary conditions. At eta = 0,
% a reflecting boundary is fine b/c mu_eta > 0 at eta = 0
% and sigma_eta -> 0. At eta = etastar, the reflecting boundary
% is an equilibrium outcome.
%
% Recall ghost nodes approach:
% (v(N+1) - v(N)) / dx = 0 => v(N+1) = v(N)
% => (v(N+1) - 2 v(N) + v(N-1)) / dx^2 = (-v(N) + v(N-1)) / dx^2;
dX = diff(eta);
N  = s.N_welfare;

L2coefs          = zeros(N, 1);
L2coefs(2:N - 1) = eta_sigma_etasq(2:N - 1) ./ (2 .* dX(1:N - 2) .* dX(2:N - 1));
L2coefs(end)     = eta_sigma_etasq(end) ./ (2 .* dX(end - 1) .* dX(end)); % result of reflecting boundary & ghost nodes approach
if (isfield(s, 'tail_risk') && s.tail_risk) || (isfield(s, 'issuance') && s.issuance > 1)
    L2coefs(1)   = eta_sigma_etasq(end) ./ (2 .* dX(end - 1) .* dX(end)); % exogenous boundary at left => do not zero vol out
end

DU = zeros(N, 1);
DD = zeros(N, 1);
D0 = zeros(N, 1);

% DU(2:N) line: Upwind fd => use previous node's drift + vol term
% DD(1:N-1) line: Upwind fd => use next node's drift + vol term
% D0(1:N-1) line: want forward diff => subtract current node value, adjusted for drift and vol term at current node.
%                 Note that this also subtracts the vol term. Finally,
%                 b/c we do not fill in D0(N) with a term from DU,
%                 if DD(N-1) was just the vol term (no drift b/c mu_eta > 0),
%                 then DD(N-1) = vol term and D0(N) = vol term
%                 => vol term * (V(N - 1) - V(N)) at the boundary, as desired by reflecting boundary
%                 e.g. Poisson equation: -u''(x) = 0 w/ u'(x_{end}) = 0 would result
%                      in a stencil of [1, -1] at indices (N, N-1) and (N, N), resp. b/c
%                      u''(x) \approx (u(N+1) - 2 u(N) + u(N-1)) / 2 dx^2 = (-u(N) + u(N-1)) / 2dx^2
% D0(2:N) line: backward diff => "add" current node value, adjusted for drift and vol term
%               Note that this also subtracts the vol term once more for a total of 2 subtractions from D0
%               (and recall centered second difference is (v(N+1) - 2 v(N) + v(N-1) / 2 dX^2))
% Last line: account for the discount rate term from HJB
DU(2:N)     = max(eta_mu_eta(1:N - 1), 0) ./ dX + L2coefs(1:N - 1);
DD(1:N - 1) = max(-eta_mu_eta(2:N), 0.) ./ dX + L2coefs(2:N);
D0(1:N - 1) = -DU(2:N);
D0(2:N)     = D0(2:N) - DD(1:N - 1);
D0          = D0 - s.r .* ones(N, 1);

A = spdiags([DD D0 DU], -1:1, N, N);

% To check
% DU1 = zeros(N, 1); DU2 = zeros(N, 1);
% DD1 = zeros(N, 1); DD2 = zeros(N, 1);
% D01 = zeros(N, 1); D02 = zeros(N, 1);
% DU1(2:N) = max(eta_mu_eta(1:N - 1), 0) ./ dX;
% DU2(2:N) = L2coefs(1:N - 1);
% DD1(1:N - 1) = max(-eta_mu_eta(2:N), 0.) ./ dX;
% DD2(1:N - 1) = L2coefs(2:N);
% D01(1:N - 1) = -DU1(2:N);
% D02(1:N - 1) = -DU2(2:N);
% D01(2:N) = D01(2:N) - DD1(1:N - 1);
% D02(2:N) = D02(2:N) - DD2(1:N - 1);
% L1 = spdiags([DD1, D01, DU1], -1:1, N, N);
% L2 = spdiags([DD2, D02, DU2], -1:1, N, N);
% Acheck = L1 + L2 - spdiags(s.r * ones(N, 1), 0, N, N);

% Inhomogeneous parts
muK        = investment_fnct(Q, 'Phi', s) - s.delta;
conv_yield = convenience_yield(x_hd, 0, s);
pi_costs   =  s.lambda / 2 * (inflation - s.pi_target).^2;
EdlogK     = (muK - s.sigma^2 / 2);
Inhom      = log(s.r * Q .* (1 + (theta - 1) .* eta)) ...
    + conv_yield - pi_costs + EdlogK / s.r;

% Solve linear equation system to get welfare
% e.g. for total welfare,
% ((eta_sigma_eta)^2 / 2) * L2 * H + eta_mu_eta * L1 * H - r * H = -Inhom
unif_H   = A \ (-Inhom);
unif_Vd  = A \ (-conv_yield);
unif_Vpi = A \ (-pi_costs);

% Interpolate back to non-uniform grid
% and note that
% V = log(Q(1 + (theta - 1) * eta)) / r + h = Vc + Vd - Vpi
H   = interp1(eta, unif_H,   orig_eta);
Vd  = interp1(eta, unif_Vd,  orig_eta);
Vpi = interp1(eta, unif_Vpi, orig_eta);
V   = H;
Vc  = V - (Vd - Vpi);

end
