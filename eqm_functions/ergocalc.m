function [etaout, densityout, CDF, enudge] = ergocalc(eta_in, Dyn, s)
% [etaout, densityout, CDF, enudge] = ergocalc(eta_in, Dyn, s)
%
% This function calculates the ergodic (stationary) distribution.
%
% Written by William Chen and Gregory Phelan, Jan. 2022

eta_mu_eta    = Dyn(:, 4) .* eta_in;
eta_sigma_eta = Dyn(:, 2) .* eta_in;
odefun = @(eta,f) ergofnct(eta, f, eta_in, eta_mu_eta, eta_sigma_eta, s);
options = odeset('events', 'evntergo', 'NonNegative', 1);

% we solve starting from \eta^*
[e1, f1] = ode45(odefun, flipud(eta_in), 100, options); % Any initial value works, normalize later

d1 = f1 ./ (interp1(eta_in, eta_sigma_eta, e1, s.interp).^2);

e1 = e1((d1 < Inf));
d1 = d1((d1 < Inf));

if isfield(s, 'smooth_stationary_pdf')
    if s.smooth_stationary_pdf
        smooth_i = 1:s.smooth_stationary_pdf_N:length(e1);
        d1 = interp1(e1(smooth_i), d1(smooth_i), e1, s.smooth_stationary_pdf_interp);
    end
end

% integrate to 1, Riemann sum approximation is fine
int = abs(d1(2:end)' * diff(e1));

d = d1(2:end) ./ int;

%get CDF
ee = flip(e1);
dd = flip(d);
etaout = ee(1:end - 1);
densityout = dd;

% diffee = diff(ee);
% TT  = tril(ones(length(diffee)));
% CDF = TT * (dd .* diffee);
pdf_interp = griddedInterpolant(etaout, densityout);
CDF = etaout;
for i = 1:length(etaout)
    CDF(i) = integral(@(x) pdf_interp(x), etaout(1), etaout(i));
end

if s.massdiff ~= 0
    massdiff=s.massdiff;
    cdffun=@(eta) interp1(ee(2:end),CDF,eta,s.interp);
    dataset(cdffun(s.center),cdffun(s.gstart))
    if cdffun(s.center)-cdffun(s.gstart)>massdiff
    enudge=fzero(@(nud) cdffun(s.center)-cdffun(s.center-nud)-massdiff,[0 .99*s.center]);
    else
        enudge=s.center-s.gstart;
        disp('massdiff too big!')
    end
else
    enudge=0;
end
end
