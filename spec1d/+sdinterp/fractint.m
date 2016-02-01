function [xi,yi,ci] = fractint(x, y, xp, s, tol)

% FRACTINT Self-Affine Fractal Interpolation Function
%    [XI,YI,CI] = FRACTINT(X,Y,XP,S,TOL) employs an affine map within
%    an Iterated Function System (IFS) to interpolate YI, the values
%    of the underlying function Y, at the points in the array XI.  XP
%    specifies the locations of the interpolation points.  X and Y
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
%    TOL specifies the relative error of the SA-FIF convergence.
%    It defaults to 1E-4.
%
%    Example:
%       x = [-1 0 1 2 3];
%       y = [-1 0 1 8 9];
%       s = [-0.5 0.23 -0.51080 -0.35670];
%       xp = 0:0.1:1;
%       [xi,yi,ci] = fractint(x, y, xp, s);
%       plot(x,y,'ko',xi,yi,'r.')
%       legend('Original data','Interpolated data')

% Joe Henning - Fall 2013

% B.F. Barnsley, J. Elton, D. Hardin, and P. Massopust
% Hidden variable fractal interpolation functions
% SIAM J. Math. Ana. 20:1218-1242
% 1989

n = length(x);

if (n ~= length(y))
   fprintf('??? Bad y input to fractint ==> y must be length(x)\n');
   xi = [];
   yi = [];
   ci = [];
   return
end

if (nargin < 4)
   s = zeros(1,n-1);
   tol = 1e-4;
elseif (nargin < 5)
   tol = 1e-4;
end

if ((n-1) ~= length(s))
   fprintf('??? Bad s input to fractint ==> s must be length(x)-1\n');
   xi = [];
   yi = [];
   ci = [];
   return
end

for i = 1:length(s)
   if (abs(s(i)) >= 1)
      fprintf('??? Bad s input to fractint ==> abs(s) < 1\n');
      xi = [];
      yi = [];
      ci = [];
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

xi = [];
yi = [];
ci = [];
for i = 1:length(xp)
   newx = -999;
   xtemp = x(2);
   ytemp = y(2);
   cnt = 0;
   if (abs(xp(i)) < tol)
      pfrac = 1;
   else
      pfrac = xp(i);
   end
   while (abs((newx-xp(i))/pfrac) > tol)
      cnt = cnt + 1;
      m = fix(2 + (n-1)*rand);
      newx = a(m)*xtemp + b(m);
      newy = c(m)*xtemp + s(m)*ytemp + d(m);
      xtemp = newx;
      ytemp = newy;
   end
   ci = [ci; cnt];
   xi = [xi; xtemp];
   yi = [yi; ytemp];
end
