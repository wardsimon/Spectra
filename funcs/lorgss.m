function [y, name, pnames, pin]=lorgss(x,p, flag)
% lorgss    : Lorentzian + gaussian
% [y, {name, pnames, pin}]=lorgss(x,p, {flag})
%
% MFIT Lorentzian + gaussian fitting function
% p = [ L_amp L_centre L_width Gauss_amp Gauss_centre Gauss_width Background ]

% Author: MZ <mzinkin@sghms.ac.uk> 
% Description:  Lorentzian + gaussian

if nargin==2
	a1=p(1);
	a2=p(4);
	c1=p(2);
	c2=p(5)+p(2);
	w1=p(3);
	w2=p(6);
	bg=p(7);

	y=p(7) + p(1)./(1+ (x-p(2)).^2/p(3)^2) + ...
	     p(4)*exp(-0.5*((x-p(5))/p(6)).^2) ;
else
	y=[];
	name='Lorentzian + Gaussian';
	pnames=str2mat('L.amplitude','L.centre','L.width',...
	              'Gauss.amplitude','Gauss.centre',...
	    'Gauss.width','Background');
	if flag==1, pin=[0 0 1 0 0 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on lorentzian peak');
		[pin(2) pin(1)]=ginput(1);
		mf_msg('Click on lorentzian width');
		[pin(3) y]=ginput(1);
		pin(3)=abs(pin(3)-pin(2));
		mf_msg('Click on gaussian peak');
		[pin(5) pin(4)]=ginput(1);
		mf_msg('Click on gaussian width');
		[pin(6) y]=ginput(1);
		pin(6)=abs(pin(6)-pin(5));
		mf_msg('Click on background');
		[x pin(7)]=ginput(1);
		pin(1)=pin(1)-pin(7);
		pin(4)=pin(4)-pin(7);
	end
end

