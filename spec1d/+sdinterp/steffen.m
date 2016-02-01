function yi = steffen(x, y, xi, yp)

% STEFFEN 1-D Steffen interpolation
%    STEFFEN(X,Y,XI) interpolates to find YI, the values of the
%    underlying function Y at the points in the array XI, using
%    the method of Steffen.  X and Y must be vectors of length N.
%
%    Steffen's method is based on a third-order polynomial.  The
%    slope at each grid point is calculated in a way to guarantee
%    a monotonic behavior of the interpolating function.  The 
%    curve is smooth up to the first derivative.

% Joe Henning - Summer 2014

% M. Steffen
% A Simple Method for Monotonic Interpolation in One Dimension
% Astron. Astrophys. 239, 443-450 (1990)

if nargin < 4
   yp = [];
end

n = length(x);

if (isempty(yp))
   % calculate slopes
   yp = zeros(1,n);

   for i = 2:n-1
      hi = x(i+1) - x(i);
      him1 = x(i) - x(i-1);
      si = (y(i+1) - y(i))/hi;
      sim1 = (y(i) - y(i-1))/him1;
      pi = (sim1*hi + si*him1)/(him1 + hi);
      
      if (sim1*si <= 0)
         yp(i) = 0;
      elseif (abs(pi) > 2*abs(sim1))
         if (sign(sim1) ~= sign(si))
            error('Sign si not equal to sim1');
         end
         a = sign(sim1);
         yp(i) = 2*a*min([abs(sim1),abs(si)]);
      elseif (abs(pi) > 2*abs(si))
         if (sign(sim1) ~= sign(si))
            error('Sign si not equal to sim1');
         end
         a = sign(sim1);
         yp(i) = 2*a*min([abs(sim1),abs(si)]);
      else
         yp(i) = pi;
      end
   end

   % first point
   h1 = x(2) - x(1);
   h2 = x(3) - x(2);
   s1 = (y(2) - y(1))/h1;
   s2 = (y(3) - y(2))/h2;
   p1 = s1*(1 + h1/(h1 + h2)) - s2*h1/(h1 + h2);
   if (p1*s1 <= 0)
      yp(1) = 0;
   elseif (abs(p1) > 2*abs(s1))
      yp(1) = 2*s1;
   else
      yp(1) = p1;
   end

   % last point
   hnm1 = x(n) - x(n-1);
   hnm2 = x(n-1) - x(n-2);
   snm1 = (y(n) - y(n-1))/hnm1;
   snm2 = (y(n-1) - y(n-2))/hnm2;
   pn = snm1*(1 + hnm1/(hnm1 + hnm2)) - snm2*hnm1/(hnm1 + hnm2);
   if (pn*snm1 <= 0)
      yp(n) = 0;
   elseif (abs(pn) > 2*abs(snm1))
      yp(n) = 2*snm1;
   else
      yp(n) = pn;
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
      fprintf('??? Bad x input to steffen ==> x values must be distinct\n');
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

   s = (y(khi) - y(klo))/h;

   a = (yp(klo) + yp(khi) - 2*s)/h/h;
   b = (3*s - 2*yp(klo) - yp(khi))/h;
   c = yp(klo);
   d = y(klo);

   t = xi(i) - x(klo);
   t2 = t*t;
   t3 = t2*t;

   yi(i) = a*t3 + b*t2 + c*t + d;
end
