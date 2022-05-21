function fp=ergofnct(eta,f,etaout,mu_eta,sigma_eta,s)

mu=interp1(etaout,mu_eta,eta,s.interp);
sigma=interp1(etaout,sigma_eta,eta,s.interp);

fp=2*mu*f/(sigma^2);
