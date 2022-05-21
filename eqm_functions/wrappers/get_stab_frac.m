function [out1, out2, out3, out4] = get_stab_frac(input, s, input_type)
% This function takes in parameters governing the desired monetary policy
% rule and returns the stability fraction.
%
% Written by William Chen and Greg Phelan, Aug. 2020

if strcmp(input_type, 'mp')
    if s.nomruletype == 1
        if length(input) > 2
            s.iput = input(1);
            s.ireg = input(2);
            s.ilaw = input(3);
        else
            s.iput = input(1);
            s.ilaw = input(2);
        end
    elseif s.nomruletype == 0
        s.iput = input(1);
        s.ilaw = input(2);
    end
elseif strcmp(input_type, 'eta')
    if length(input) == 1
        s.etaLAW = input(1);
    elseif length(input) == 2
        s.etaPUT = input(1);
        s.etaLAW = input(2);
    end
else
    error('The input_type must be either mp or eta.');
end

[sol, summary_stats] = geteqm(s);

out1 = summary_stats.ergodic.stab_frac;
out2 = sol.eta(find(sol.psi == 1, 1) - 1);

if nargout > 2
    out3 = sol.interestvec;
    out4 = sol.eta;
end
end
