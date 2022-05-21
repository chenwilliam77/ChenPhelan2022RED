function [value, isterminal, direction] = evntergo(eta, F)
% [value,isterminal,direction] = evntergo(eta, F) returns three
% variables used by ode45 to determine when to terminate integration when
% calculating ergodic distribution


value = [F];  %
isterminal = [1 ];   %
direction = [-1 ];    % event occurs when we get there from above
