function [y, name, pnames, pin]=dirac(x,p, flag)
% dirac     : Dirac peak
% [y, {name, pnames, pin}]=Dirac(x,p, {flag}) 
% MFIT Dirac fitting function
% p = [ Amp Position ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
    y = 0*x;
    i1 = find(x<= p(2));  if isempty(i1), i1 = length(x); end
    i1 = i1(end);
    i2 = find(x > p(2));  if isempty(i2), i2 = 1; end
    i2 = i2(1);
    if abs(x(i1) - p(2)) > abs(x(i2) - p(2))
	y(i2) = p(1);
    else
	y(i1) = p(1);
    end
else
	y=[];
	name='dirac';
	pnames=str2mat('Amplitude','position');
	if flag==1, pin=[1 1]; else pin = p; end
	if flag==2
		mf_msg('Click on position');
		[cen amp]=ginput(1);

		pin=[ amp cen ];
	end
end
