function [y, name, pnames, pin] = airyfp(x,p, flag)
% airyfp    : Airy function for Fabry Perot
% [y, name, pnames, pin] = airyfp(x,p, flag)
%  = I/(1+M*(sin(Phi)^2))+bg  with Phi = pi*(x-x0)/Period
%
% Airy function for Fabry-Perot with
%    p = [ x0 M Intensity Period bg]

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Airy function for Fabry-Perot

if nargin==2
	y = (p(3) ./ (1 + p(2) * (sin( pi*(x-p(1))/p(4) ) ).^2)) + p(5);
else
	y=[];
	name='Airy function';
	pnames=str2mat('Centre','M(width)', 'Intensity', 'Period(FSR)','Background');
	if flag==1, pin=[0 10 1 100 0]; else pin = p; end
	if flag==2
		mf_msg('Click on background');
		[cen bg]=ginput(1);
		mf_msg('Click on central peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on second peak for Free Spectral Range.');
		[FSR y]=ginput(1);
		FSR = abs(FSR - cen);
		pin=[ cen 1/width amp FSR bg];
	end
end
