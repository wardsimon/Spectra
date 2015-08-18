function [y, name, pnames, pin]=lorzn(x,p, flag)
% lorzn     : Lorentzian Power n
% [y, {name, pnames, pin}]=lorzn(x,p, {flag})
%
% MFIT Lorentzian power n fitting function
% p = [ Amplitude Centre Width Power Background ]

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Lorentzian Power n

if nargin==2
    y=p(5)+p(1) ./ ((1+ (x-p(2)).^2/p(3)^2 ).^p(4));
else
	y=[];
	name='Lorentzian Power N';
	pnames=str2mat('Amplitude','Centre','Width','Power','Background');
	if flag==1, pin=[0 0 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width 1 bg];
	end
end

