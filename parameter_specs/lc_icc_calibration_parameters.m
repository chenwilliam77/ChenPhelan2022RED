% This script defines the parameters used in the model.
% This script also specifies an ICC leverage constraint policy

calibration_parameters;

% Macropru parameters
lc.type = 'icc';
lc.alpha = 2;
s.leverage_constraint = lc;
