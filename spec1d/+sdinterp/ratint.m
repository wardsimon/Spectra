function [yi, R] = ratint(x, y, xi)

% RATINT Interpolation using simple rational polynomials
%    RATINT(X,Y,XI) interpolates to find YI, the value of
%    the underlying function Y at the points in the array XI,
%    using rational polynomials.  X and Y must be vectors of
%    length N.
%
%    [YI,R] = RATINT() also returns the rational polynomial
%    table R calculated for the last XI.

% Joe Henning - Fall 2011

n = length(x);

tol = n * max(y) * eps(class(y));

for k = 1:length(xi)
   xd = [];
   for i = 1:n
      xd(i) = abs(x(i) - xi(k));
   end

   [xds,i] = sort(xd);

   x = x(i);
   y = y(i);

   R = zeros(n,n);
   R(:,1) = y(:);

   % Evaluate rational polynomial
   for i = 1:n-1
      for j = 1:(n-i)
         R(j,i+1) = (x(j+i)-x(j))/((xi(k)-x(j))/(R(j+1,i)+tol) + (x(j+i)-xi(k))/(R(j,i)+tol));
      end
   end

   yi(k) = R(1,n);
end
