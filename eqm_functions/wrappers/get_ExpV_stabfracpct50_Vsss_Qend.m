function [ExpV, stab_frac, Vsss, Qend] = get_ExpV_stabfracpct50_Vsss_Qend(input, s, input_type)
% This function takes in parameters governing the desired monetary policy
% rule and returns the expected welfare and stability fraction.
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
elseif strcmp(input_type, 'constpi')
    [s, converged] = find_constant_inflation_mp_rule(s);
    if ~converged
        error('Convergence failed: constant inflation rule not found.');
    end
else
    error('The input_type must be mp, eta, pw, fp, or constpi.');
end

[sol, summary_stats] = geteqm(s);

ExpV = summary_stats.ergodic.ExpV;
stab_frac = summary_stats.ergodic.stab_frac_pct50;
Vsss = summary_stats.sss.Welfare(end); % by default, V(final eta*)
Qend = sol.Q(end);
end
