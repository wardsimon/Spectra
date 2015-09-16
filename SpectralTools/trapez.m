function q = trapez(x,y,a,b)
% trapez : Trapezoidal integration
%Syntax: q = trapez(x,y,{a,b})   or q = trapez(y)
%
% Computes the integral value of y(x) between limits a and b by trapezoidal rule.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  trapezoidal integration

% Written by: Duane Hanselman, University of Maine, (207)-581-2246
% revised by E.Farhi.

if (nargin <1)
	error('q = trapez(x,y,{a,b})   or q = trapez(y)')
end

if (nargin == 1)
	y = x;
	x = 1:length(y);
end

x = x(:); y = y(:);   % make sure x and y are row vectors

ly = length(y);
lx = length(x);
if (lx ~= ly)
  error(' x and y must be the same length')
end

if (nargin <=2 )
	a=x(1);
	b=x(lx);
end

q=find( (x>=a) & (x<=b) );
x=x(q);
y=y(q);
ly = length(y);

q = ( [0 y] + [y 0] ) / 2 .* ( [x 0] - [0 x] );
q = sum( q(2:ly) );
