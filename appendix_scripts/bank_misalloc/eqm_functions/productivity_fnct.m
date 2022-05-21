function out = productivity_fnct(eta, order, s)
% out = productivity_fnct(eta, order, s)
%
% Function determining the productivity of banks
% specified as a function of eta.
%
% eta: state
% order: 0 for level, 1 for derivative
% s: parameters and settings struct
%
% Written by William Chen and Greg Phelan, Mar. 2022


if isfield(s, 'bank_misalloc')
    bm = s.bank_misalloc; % struct with parameters related to bank misallocation
else
    bm.eta_min = 0;
    bm.eta_max = 0;
    bm.ab_min = s.ab;
end

 % Linearly decreasing rule
if order == 0
	if (eta > bm.eta_min) && (eta < bm.eta_max)
	    out = s.ab - (s.ab - bm.ab_min) * (eta - bm.eta_min) / (bm.eta_max - bm.eta_min);
    elseif eta >= bm.eta_max
        out = bm.ab_min;
	else
	    out = s.ab;
	end
elseif order == 1
	if (eta > bm.eta_min) && (eta < bm.eta_max)
	    out = - (s.ab - bm.ab_min) / (bm.eta_max - bm.eta_min);
	else
	    out = 0;
	end
else
    error('The input argument order must be 0 or 1');
end

end
