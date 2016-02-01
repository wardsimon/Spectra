function [xi,yi] = safif(x, y, s, k)

% SAFIF Self-Affine Fractal Interpolation Function
%    [XI,YI] = SAFIF(X,Y,S,K) employs an affine map within an Iterated
%    Function System (IFS) to interpolate YI, the values of the
%    underlying function Y, at the points in the array XI.  X and Y
%    must be vectors of length N.
%
%    IFS:
%        [x]   [ai  0][x]   [bi]
%     wi [y] = [ci si][y] + [di]
%
%    S specify the vertical scaling parameters associated with the
%    the mapping.  S must be a vector of length N-1.  It defaults to
%    a vector of zeros.
%
%    The S parameters must also satisfy
%        |Si| < 1 for i = 1, 2, ..., N-1
%    such that the graph of the interpolating function interpolate Y.
%    The interpolating function is called the Self-affine Fractal
%    Interpolation Function (SA-FIF).
%
%    K specifies the number of interpolation points used to form
%    XI and YI.  It defaults to 10000.
%
%    Example:
%       x = [-1 0 1 2 3];
%       y = [-1 0 1 8 9];
%       s = [-0.5 0.23 -0.51080 -0.35670];
%       [xi,yi] = safif(x, y, s);
%       plot(x,y,'ko',xi,yi,'r.')
%       legend('Original data','Interpolated data')

% Joe Henning - Fall 2013

% B.F. Barnsley, J. Elton, D. Hardin, and P. Massopust
% Hidden variable fractal interpolation functions
% SIAM J. Math. Ana. 20:1218-1242
% 1989

n = length(x);

if (n ~= length(y))
   fprintf('??? Bad y input to safif ==> y must be length(x)\n');
   xi = [];
   yi = [];
   return
end

if (nargin < 3)
   s = zeros(1,n-1);
   k = 10000;
elseif (nargin < 4)
   k = 10000;
end

if ((n-1) ~= length(s))
   fprintf('??? Bad s input to safif ==> s must be length(x)-1\n');
   xi = [];
   yi = [];
   return
end

for i = 1:length(s)
   if (abs(s(i)) >= 1)
      fprintf('??? Bad s input to safif ==> abs(s) < 1\n');
      xi = [];
      yi = [];
      return
   end
end

delx = x(n) - x(1);
s = [0 s];

a = [0];
b = [0];
c = [0];
d = [0];
for i = 2:n
   a(i) = (x(i) - x(i-1))/delx;
   b(i) = (x(n)*x(i-1) - x(1)*x(i))/delx;
   c(i) = (y(i) - y(i-1) - s(i)*(y(n) - y(1)))/delx;
   d(i) = (x(n)*y(i-1) - x(1)*y(i) - s(i)*(x(n)*y(1) - x(1)*y(n)))/delx;
end

xtemp = x(2);
ytemp = y(2);
xi = [xtemp];
yi = [ytemp];
for i = 1:k
   m = fix(2 + (n-1)*rand);
   newx = a(m)*xtemp + b(m);
   newy = c(m)*xtemp + s(m)*ytemp + d(m);
   xtemp = newx;
   ytemp = newy;
   xi = [xi; newx];
   yi = [yi; newy];
end

% box-counting dimension of the interpolating function
% for equidistant points and sum(abs(s)) > 1
%boxdim = 1 + log(sum(abs(s)))/log(n-1)
