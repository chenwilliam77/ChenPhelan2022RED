function out = convenience_yield(x_hd, order, s)
% out = convenience_yield(x_hd, order, s)
%
% Implements the convenience yield function
%
% v(x_hd) = beta_1 * (x_hd + beta_2) ^ (1 - beta_3) / (1 - beta_3)
%
% x_hd  - share of deposits to households' non-bank net worth,
% order - 0 for level, 1 for first derivative,
%         results in an error otherwise
%
% Written by William Chen and Gregory Phelan, Jan. 2022

if order == 0
    if s.beta_3 == 1
		out = s.beta_1 .* log(x_hd + s.beta_2);
    else
		out = s.beta_1 .* (x_hd + s.beta_2) .^ (1 - s.beta_3) ./ (1 - s.beta_3);
    end
elseif order == 1
    out = s.beta_1 .* (x_hd + s.beta_2) .^ (-s.beta_3);
else
	error("The second argument `order` must be a 0 or a 1");
end


end
