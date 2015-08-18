function [y, name, pnames, pin]=polynomial(x,p, flag)
% polynomial: polynomial
% [y, {name, pnames, pin}]=polynomial(x,p, {flag})
%
% MFIT polynomial function 
%  y = p(1)*x^d + p(2)*x^(d-1) + ... + p(d)*x + p(d+1)

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  polynomial

np = length(p);

if nargin==2 & np > 0;
	y=polyval(p,x);
else
	name = 'Polynomial';
	pin = p;
	y=zeros(size(x));
	if flag==2
		np = 0;	
		but = 1;
		x=[]; y=[];
		pin = [];
		while but==1
			mf_msg(sprintf('Click on point %d (right button to end)',np+1));
			[x1 y1 but]=ginput(1);
			if but==1
				x = [ x x1 ];
				y = [ y y1 ];
				np=np+1;
			end
		end
		[ny,pin] = mf_interp(x,y,x(ceil(length(x)/2)),length(x),np-1);
	end
	y=[];
	pnames=str2mat('p1');
	if np>0 & any(p)
		name=sprintf('%dth order polynome',np);
	else
		name='n order polynome : clik on Guess button to set n.';
		if flag==1, pin=[ 0 ]; else pin = p; end
	end

	for i=2:np
		pnames=str2mat(pnames,...
		               sprintf('p%d',i));
	end
end
