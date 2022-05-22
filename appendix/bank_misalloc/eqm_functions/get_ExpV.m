function out = get_ExpV(input, s, input_type)
% This function takes in parameters governing the desired monetary policy
% rule and returns the expected welfare
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
        s.etaLAW = input(1) + input(2);
    end
elseif strcmp(input_type, 'pw')
    s.nomruletype = 1;
    if length(input) == 3
        s.iput = input(1);
        s.ilaw = input(2);
        s.etaPUT = input(3);
        s.etaLAW = s.etaPUT;
    elseif length(input) == 5
        s.iput = input(1);
        s.ireg = input(2);
        s.ilaw = input(3);
        s.etaPUT = input(4);
        s.etaLAW = input(4) + input(5);
    end
elseif strcmp(input_type, 'fp')
    s.nomruletype = 0;
    s.iput = input(1);
    s.ilaw = input(2);
    s.etaPUT = input(3);
    s.etaLAW = input(3) + input(4);
else
    error('The input_type must be mp, eta, pw, or fp.');
end

[~, summary_stats] = geteqm(s);

out = summary_stats.ergodic.ExpV;
end
