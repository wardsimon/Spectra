function [ y, name, pnames, pin ] = couplph(x, p, flag)
% couplph   : Phonon coupled with pseudo spin system.
% [ y, name, pnames, pin ] = couplph(x, p, flag)
% This function describes a phonon mode coupled with a
% dynamical pseudo-spin system. 
% It should specially suit non centro-symetric ferroelectric crystals.
% Parameters are :
% p(1) : omega phonon (Wo)
% p(2) : Amplitude
% p(3) : Temperature (T)
% p(4) : Transition temperature (Tc)
% p(5) : omega spin flip
% p(6) : tau tilde = spin-phonon coupling (0=strong)
% p(7) : background



if nargin == 2
	w0 = p(1);
	T = p(3);
	Tc = p(4);
	wf = p(5); wt = w0/wf;
	tt = p(6);

	tau = (T-Tc)/Tc;
	w = x/w0;
	tmp = tau+1-tt;

	% spin-spin correlation function (spectra)
	% yss = wt/w0*(tau+1)/sqr(tmp) * sqr(sqr(w)-1) ./ ( sqr(sqr(w)-tau/tmp) + sqr(wt*w.*(sqr(w)-1)/tmp) );

	% phonon-phonon correlation function (spectra)
	yqq = Tc*(1-tt)/sqr(w0) * wt/w0*(tau+1)/sqr(tmp) ./ ( sqr(sqr(w)-tau/tmp) + sqr(wt*w.*(sqr(w)-1)/tmp) );

	y = p(2)*yqq + p(7);
else
	name = 'Pseudospin-phonon coupled';
	y = [];
	pnames = str2mat(...
		'Phonon energy (Wo)',...
		'Amplitude',...
		'Temperature (T)',...
		'Transition temp. (Tc)',...
		'Spin energy (Wf)',...
		'Tau (coupling)',...
		'Background');
	if flag == 1
		pin = [  1 0 25 30 1 0.5 1 ];
	else
		pin = p;
	end
	if flag == 2
		mf_msg('Click on peak');
		[w0 amp]=ginput(1);
		pin(1) = w0;
		pin(2) = amp/p(4);
	end
end
