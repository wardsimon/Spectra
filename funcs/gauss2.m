function [y, name, pnames, pin]=gauss2(x,p, flag)
% gauss2    : Gaussian Squared
% [y, {name, pnames, pin}]=gauss2(x,p, {flag}) 
%
% MFIT Gaussian Squared fitting function
% p = [ Amp Centre Width BachG ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Gaussian Squared

if nargin==2;
    y=p(4)+p(1)*exp(-0.5*((x-p(2))/p(3)).^2).^2;
else
	y=[];
	name='Gaussian Squared';
	pnames=str2mat('Amplitude','Centre','Width','Background');
	if flag==1, pin=[0 0 1 1]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width bg];
	end
end
