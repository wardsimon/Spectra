function [y, name, pnames, pin] = bose(x,p, flag)
% bose      : bose factor
% y = bose(x,p)
%   = 1/(exp(p*x) - 1)
%
% Bose function with
%    p = h/2pi/Kb/T   in 'x' units

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Bose function


if nargin==2
	y = 1 ./ (exp(p(1) * x) - 1);
else
	y=[];
	name='bose';
	pnames=str2mat('Tau(h/kT)');
	if flag==1, pin=[1]; else pin = p; end
	if flag==2
		mf_msg('No guess');
		pin=[ 1 ];
	end

end
