% This script defines the parameters used in the model.
% This script also specifies an ICC leverage constraint policy

calibration_parameters;

% Moral Hazard parameters
bm.eta_min = 0.05;
bm.eta_max = 0.05;
bm.ab_min = s.ab * 0.95;
s.bank_misalloc = bm;
