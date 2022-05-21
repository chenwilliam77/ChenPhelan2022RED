function out = interest(eta, s)
% This function outputs the central bank's interest targeting rule.
% The rule is subject to a non-negativity constraint.
%
% Written by William Chen and Greg Phelan, Aug 2020

if s.nomrule == 1 % Piecewise or fed put style rule
    imax = s.imax;

    if s.etaPUT > s.etaLAW
        error('etaPUT cannot be larger than etaLAW');
    end

    if s.nomruletype == 1
        % piecewise rule--add two breakpoints
        out = min(s.ireg, imax);

        if (eta < s.etaPUT)
            out = s.iput;
        end
        if (eta > s.etaLAW)
            out = min(s.ilaw, imax);
        end

    else
        % Fed Put, from iLAW to iPUT with linear rule in between
        ireg = (s.ilaw - s.iput) / (s.etaLAW - s.etaPUT) * (eta - s.etaPUT) + s.iput;

        out = min(ireg, imax);

        if (eta < s.etaPUT)
            out = s.iput;
        end
        if (eta > s.etaLAW)
            out = min(s.ilaw, imax);
        end
    end
elseif s.nomrule == 2
    % Function specified by the field nomrule_mp_interp
    % in the settings struct
    out = s.nomrule_mp_interp(eta);
else
    out = 0;
end
end
