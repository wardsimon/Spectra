function [y, name, pnames, pin]=lorz2(x,p, flag)
% lorz2     : Lorentzian squared
% [y, {name, pnames, pin}]=lorz2(x,p, {flag})
%
% MFIT Lorentzian squared fitting function
% p = [ Amplitude Centre Width Background ]


% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Lorentzian squared

if nargin==2
    y=p(4)+p(1) ./ ((1+ (x-p(2)).^2/p(3)^2 ).^2);
else
	y=[];
	name='Lorentzian Squared';
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

