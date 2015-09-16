function [y, name, pnames, pin]=triangl(x,p, flag)
% triangl   : Triangular
% function [y, {name, pnames, pin}]=triangl(x,p, {flag})
%
% MFIT Triangular fitting function
% p = [ Amplitude Centre Width Background ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Triangular

if nargin==2;
	s=sign(x-p(2));
	y=(p(3)-s.*(x-p(2)))/p(3)^2;
	f=find(y<0);
	y(f)=zeros(size(f));
	y=y*p(1)+p(4);
else
	y=[];
	name='Triangle';
	pnames=str2mat('Amplitude','Centre','Width','Background');
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
		pin=[amp*width cen width bg];
	end
end
