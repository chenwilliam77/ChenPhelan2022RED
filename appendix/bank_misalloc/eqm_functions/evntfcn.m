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
