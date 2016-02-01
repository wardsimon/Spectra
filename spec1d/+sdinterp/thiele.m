function [yi, T] = thiele(x, y, xi, p)

% THIELE Interpolation using rational polynomials
%    THIELE(X,Y,XI,P) interpolates to find YI, the values of
%    the underlying function Y at the points in the array XI,
%    using rational polynomials.  X and Y must be vectors of
%    length N.
%
%    P specifies the degree of the numerator.  It defaults
%    to a value of 1 (Thiele interpolation).
%
%    [YI,T] = THIELE() also returns the rational polynomial
%    table T calculated for the last XI.

% Some Techniques for Rational Interpolation
% F. M. Larkin
% The Computer Journal (1967), 10 (2)
% pgs. 178-187

% Joe Henning - Fall 2014

if (nargin < 4)
   p = 1;
end

n = length(x);

if (p < 1 || p >= (n-1))
   fprintf('??? Bad p input to thiele ==> 1 <= p < (n-1)\n');
   yi = NaN;
   T = [];
   return
end

q = (n - 1) - p;

tol = n * max(y) * eps(class(y));

for k = 1:length(xi)
   xd = [];
   for i = 1:n
      xd(i) = abs(x(i) - xi(k));
   end

   [xds,i] = sort(xd);

   x = x(i);
   y = y(i);

   T = zeros(n,n);
   T(:,1) = y(:);

   % Evaluate rational polynomial
   for i = 1:n-1
      for j = 1:(n-i)
         if (i <= p)
            T(j,i+1) = ((xi(k)-x(j))*T(j+1,i) + (x(j+i)-xi(k))*T(j,i))/(x(j+i)-x(j));
         else
            T(j,i+1) = T(j+1,i-1) + (x(j+i)-x(j))/((xi(k)-x(j))/(T(j+1,i)-T(j+1,i-1)+tol) + (x(j+i)-xi(k))/(T(j,i)-T(j+1,i-1)+tol));
         end
      end
   end

   yi(k) = T(1,n);
end
