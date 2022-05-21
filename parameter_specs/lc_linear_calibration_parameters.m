% This script defines the parameters used in the model.
% This script also specifies a linear leverage constraint policy

calibration_parameters;

% Macropru parameters
lc.type = 'linear';
lc.etagrid = [0; 0.01];
lc.Ls = [Inf; Inf];
s.leverage_constraint = lc;
