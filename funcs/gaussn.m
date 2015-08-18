function [y, name, pnames, pin]=gaussn(x,p, flag)
% gaussn    : Gaussian Power n
% [y, {name, pnames, pin}]=gaussn(x,p, {flag}) 
%
% MFIT Gaussian power n fitting function
% p = [ Amp Centre Width Power BachG ]

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Gaussian Power n

if nargin==2;
    y=p(5)+p(1)*exp(-0.5*((x-p(2))/p(3)).^2).^p(4);
else
	y=[];
	name='Gaussian Power N';
	pnames=str2mat('Amplitude','Centre','Width','Power','Background');
	if flag==1, pin=[0 0 1 1 1]; else pin = p; end
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
