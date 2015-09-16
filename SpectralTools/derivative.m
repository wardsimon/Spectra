function [deriv,dy] = derivative(x,y,order,last)
% derivative : Derivative approximation
%Syntax: [deriv, dif_table] = derivative(x,y,{order=1,mode='silent'}) or derivative(y)
%
% Computes the numerical derivative of 'y' versus 'x' (sorted)
% using Gregory-Newton forward formula at order specified.
% also return finite differences table.
% NB : It is advised to smooth data before derivating.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  derivative approximation by Gregory-Newton differences formula.

% uses :
% Part of 'Spectral tools'. E.Farhi. 12/95. From Numerical recipies in C.

% Thanks to David Cavaille for reminding me to tell users to smooth their
% data before derivating.

if (nargin < 1)
        error('usage: [deriv, dif_table] = derivative(x,y,{order = 1,mode=''silent''}) or derivative(y)');
end
if nargin < 3, order = []; end
if nargin < 4, last = []; end
if isempty(order)
  order = 1;
end

if (nargin == 1)
  y=x;
  x=1:length(y);
end

if isempty(last)
  last = 'silent';
end

if strcmp(last,'silent')
  verbose = 0;
else
  verbose = 1;
end

if (verbose>0)
  fprintf(1,'* Derivating : ');
end

y2 = y;

x = x(:);
y = y(:);

la=length(y);
order = max(min(order,la-1),1);
dy=zeros(order,la);
for n=1:(order)     % computes differences
  if (n==1)
    dy(1,1:(la-1)) = y(2:la) - y(1:(la-1));                     % order 1
  else
    dy(n,1:(la-n+1)) = dy(n-1,2:(la-n+2)) - dy(n-1,1:(la-n+1));
  end
  if (verbose>0)
    fprintf(1,'#');
  end
end

dx = x(2:la) - x(1:la-1);
dx(la) = dx(la-1);

deriv = dy(1,:);

for i=2:order
  deriv = deriv + ((-1)^(i-1))*dy(i,:)/i;
  if (verbose>0)
    fprintf(1,'#');
  end
end

deriv = deriv ./dx;
deriv = deriv(:);

deriv(la) = deriv(la-1);

if (size(y2,2) == 1)
  deriv = deriv';
  dy = dy';
end

if (verbose>0)
  fprintf(1,' ');
end

