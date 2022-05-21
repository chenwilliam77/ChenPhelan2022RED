function [L, dL] = leverage_constraint(eta, f, lc)
% [L, dL] = leverage_constraint(eta, f, lc)
%
% Implements the leverage constraint function L(eta, f(eta))
% depending on the settings specified in the struct lc.
% The first derivative must also be specified.
%
% Written by William Chen and Gregory Phelan, Aug. 2020

Ltype = lc.type;
if strcmp(Ltype, 'constant')
    L = lc.L;
    if nargout == 2
        dL = 0;
    end
elseif strcmp(Ltype, 'piecewise')
    % Find left end point of bin that eta falls into
    i = find(eta > lc.etagrid, 1, 'last');
    L = lc.Ls(i); % Find leverage over this bin
    if nargout == 2
        dL = 0;
    end
elseif strcmp(Ltype, 'linear')
    if isempty(find(0 == lc.etagrid, 1))
        error('The point 0 must be an element of lc.etagrid to use a piecewise linear leverage constraint');
    end
    if length(lc.etagrid) ~= length(lc.Ls)
        error('The length of etagrid should match the length of Ls in the leverage_constraint struct');
    end

    % Find left endpoint of bin that eta falls into
    i     = find(eta > lc.etagrid, 1, 'last');
    coef0 = lc.Ls(i); % leverage required at L(eta_i)
    coef1 = lc.Ls(i + 1); % leverage required at L(eta_{i + 1})
    if isinf(coef0)
        L = Inf;
    else
        L = (coef1 - coef0) * (eta - lc.etagrid(i)) / (lc.etagrid(i + 1) - lc.etagrid(i)) + coef0;
    end
    dL    = (coef1 - coef0) / (lc.etagrid(i + 1) - lc.etagrid(i));
elseif strcmp(Ltype, 'icc')
    L  = lc.alpha * f(1); % alpha * theta
    dL = lc.alpha * f(2); % alpha * theta'
elseif strcmp(Ltype, 'icc-constant')
    Licc  = lc.alpha * f(1); % alpha * theta
    Lcons = lc.L;

    % Choosing the more binding constraint
    if Licc <= Lcons
        L = Licc;
        dL = lc.alpha * f(2); % alpha * theta'
    else
        L  = Lcons;
        dL = 0;
    end
elseif strcmp(Ltype, 'icc-piecewise')
    Licc = lc.alpha * f(1); % alpha * theta
    i    = find(eta > lc.etagrid, 1, 'last');
    Lpw  = lc.Ls(i);

    % Choosing the more binding constraint
    if Licc <= Lpw
        L  = Licc;
        dL = lc.alpha * f(2); % alpha * theta'
    else
        L  = Lpw;
        dL = 0;
    end
elseif strcmp(lc, 'icc-linear')
    Licc  = lc.alpha * f(1); % alpha * theta
    i     = find(eta > lc.etagrid, 1, 'last');
    coef0 = lc.Ls(i); % leverage required at L(eta_i)
    coef1 = lc.Ls(i + 1); % leverage required at L(eta_{i + 1})
    Llin  = (coef1 - coef0) * (eta - lc.etagrid(i)) / (lc.etagrid(i + 1) - lc.etagrid(i));

    % Choosing the more binding constraint
    if Licc <= Llin
        L  = Licc;
        dL = lc.alpha * f(2); % alpha * theta'
    else
        L  = Llin;
        dL = (coef1 - coef0) / (lc.etagrid(i + 1) - lc.etagrid(i));
    end
elseif strcmp(Ltype, 'true-icc')
    L  = lc.alpha * f(1) - 1; % alpha * theta
    dL = lc.alpha * f(2); % alpha * theta'
elseif strcmp(Ltype, 'true-icc-constant')
    Licc  = lc.alpha * f(1) - 1; % alpha * theta
    Lcons = lc.L;

    % Choosing the more binding constraint
    if Licc <= Lcons
        L = Licc;
        dL = lc.alpha * f(2); % alpha * theta'
    else
        L  = Lcons;
        dL = 0;
    end
elseif strcmp(Ltype, 'true-icc-piecewise')
    Licc = lc.alpha * f(1) - 1; % alpha * theta
    i    = find(eta > lc.etagrid, 1, 'last');
    Lpw  = lc.Ls(i);

    % Choosing the more binding constraint
    if Licc <= Lpw
        L  = Licc;
        dL = lc.alpha * f(2); % alpha * theta'
    else
        L  = Lpw;
        dL = 0;
    end
elseif strcmp(lc, 'true-icc-linear')
    Licc  = lc.alpha * f(1) - 1; % alpha * theta
    i     = find(eta > lc.etagrid, 1, 'last');
    coef0 = lc.Ls(i); % leverage required at L(eta_i)
    coef1 = lc.Ls(i + 1); % leverage required at L(eta_{i + 1})
    Llin  = (coef1 - coef0) * (eta - lc.etagrid(i)) / (lc.etagrid(i + 1) - lc.etagrid(i));

    % Choosing the more binding constraint
    if Licc <= Llin
        L  = Licc;
        dL = lc.alpha * f(2); % alpha * theta'
    else
        L  = Llin;
        dL = (coef1 - coef0) / (lc.etagrid(i + 1) - lc.etagrid(i));
    end
else
    L = Inf;
    dL = NaN;
end

end
