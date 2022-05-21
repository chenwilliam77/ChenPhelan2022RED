function phi = consumption_equiv(V1, V0, r)
% This function calculates the consumption equivalent welfare difference
% between V1 and V0, given a discount rate r. If the output is phi,
% then giving households in the baseline scenario corresponding to V0
% (1 + phi) times more consumption would deliver the same welfare as V1.
%
% Written by William Chen and Greg Phelan, Jan. 2022
    phi = exp(r * (V1 - V0)) - 1;
end
