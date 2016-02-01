function yi = stineman(x, y, xi, yp)

% STINEMAN 1-D Stineman interpolation
%    STINEMAN(X,Y,XI) interpolates to find YI, the values of the
%    underlying function Y at the points in the array XI, using
%    the method of Stineman.  X and Y must be vectors of length N.
%
%    Stineman's method is based on a rational function which satisfies
%    three conditions outlined in his paper.  The method produces
%    monotonic curves and avoids problems seen in other interpolation
%    methods near abrupt changes of slope.  The curve is smooth up
%    to the first derivative.

% Joe Henning - Summer 2014

% Russel W. Stineman
% A Consistently Well-Behaved Method of Interpolation
% Creative Computing
% Volume 6, Number 7, 1980.
% pgs. 54-57

if nargin < 4
   yp = [];
end

n = length(x);

if (isempty(yp))
   % calculate slopes
   yp = zeros(1,n);

   for i = 2:n-1
      a = y(i) - y(i-1);
      b = y(i+1) - y(i);
      c = x(i) - x(i-1);
      d = x(i+1) - x(i);
      yp(i) = (a*(d*d + b*b) + b*(c*c + a*a))/(c*(d*d + b*b) + d*(c*c + a*a));
   end

   % first point
   s = (y(2)-y(1))/(x(2)-x(1));
   if (s > 0 & s > yp(2))
      yp(1) = 2*s - yp(2);
   elseif (s < 0 & s < yp(2))
      yp(1) = 2*s - yp(2);
   else
      yp(1) = s + abs(s)*(s - yp(2))/(abs(s) + abs(s - yp(2)));
   end

   % last point
   s = (y(n)-y(n-1))/(x(n)-x(n-1));
   if (s > 0 & s > yp(n-1))
      yp(n) = 2*s - yp(n-1);
   elseif (s < 0 & s < yp(n-1))
      yp(n) = 2*s - yp(n-1);
   else
      yp(n) = s + abs(s)*(s - yp(n-1))/(abs(s) + abs(s - yp(n-1)));
   end
end

for i = 1:length(xi)
   % Find the right place in the table by means of a bisection.
   klo = 1;
   khi = n;
   while (khi-klo > 1)
      k = fix((khi+klo)/2.0);
      if (x(k) > xi(i))
         khi = k;
      else
         klo = k;
      end
   end

   h = x(khi) - x(klo);
   if (h == 0.0)
      fprintf('??? Bad x input to stineman ==> x values must be distinct\n');
      yi(i) = NaN;
      continue;
   end
   
   isiny = 0;
   for k = 1:n
      if (xi(i) == x(k))
         yi(i) = y(k);
         isiny = 1;
         break
      end
   end

   if (isiny)
      continue
   end

   slo = (y(khi)-y(klo))/(x(khi)-x(klo));

   % ordinate corresponding to xi(i)
   y0 = y(klo) + slo*(xi(i) - x(klo));

   % vertical distances
   dylo = y(klo) + yp(klo)*(xi(i) - x(klo)) - y0;
   dyhi = y(khi) + yp(khi)*(xi(i) - x(khi)) - y0;

   % test product
   prod = dylo*dyhi;
   if (abs(prod) < eps)
      yi(i) = 0;
   elseif (prod > 0)
      yi(i) = y0 + prod/(dylo + dyhi);
   else   % prod < 0
      yi(i) = y0 + prod*(xi(i) - x(klo) + xi(i) - x(khi))/((dylo - dyhi)*(x(khi) - x(klo)));
   end
end
