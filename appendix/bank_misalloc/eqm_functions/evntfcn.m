function [value, isterminal, direction] = evntfcn(eta, F, s)
% [value,isterminal,direction] = evntfcn(eta, F) returns three
% variables used by ode45 to determine when to terminate integration.
%
% See paper for desired boundary conditions under different policies
%
% Written by William Chen and Greg Phelan, Jan. 2022

thetastar = s.thetastar;
val1 = F(1) - thetastar;

% Determine what policies apply
if isfield(s, 'dividend_subsidy')
    div_sub = s.dividend_subsidy;
else
    div_sub = false;
end

% Figure out the desired boundary conditions
if div_sub
    val2 = F(2) * eta + thetastar - 1;
else
    val2 = F(2);
end

value = [val1 val2];
isterminal = [1 1];   % terminate computation in all cases
direction = [0 0];    % event occurs whether we get there from above or below
end

% Old work on trying to implement dividend restriction
% if div_res
%    if eta < s.eta_DR % Can't issue dividends yet
%         val1 = Inf;  % Set to positive b/c theta approaches theta(etabar) from above, so event function won't terminate integration early
%         val2 = -Inf; % Set to negative number
%    elseif ~div_sub % Can issue dividends now so revert to standard boundary condition
%        if F(1) < thetastar
%            val1 = -1e-15; % val2 is going to be F(2) * eta + thetastar - 1
%        elseif F(2) * eta + thetastar - 1 > 0 % then
%            val1 = -sign(F(1)) * 3e-16;
%            val2 = -sign(F(2)) * 3e-16;
%        else
%            val2 = F(2); % normal boundary conditions hold in this case
%        end
%    end
% end

% Deprecated the Qstar and Qpstar boundary condition; these are obtained
% "for free" as a result of the theta and thetap boundary conditions
% Qstar      = s.Qmax;
% if ~exist('Qpstar_is_zero', 'var')
%     Qpstar = Qstar * (thetastar - 1) / (1 + (thetastar - 1) * eta);
% elseif Qpstar_is_zero == 0
%     Qpstar = Qstar * (thetastar - 1) / (1 + (thetastar - 1) * eta);
% elseif Qpstar_is_zero == 1
%     Qpstar = 0;
% else
% 	error('Qpstar_is_zero must be 1 or 0');
% end
% implied_Qpstar = F(3) * (F(1) - 1) / (1 + (F(1) - 1) * eta);
