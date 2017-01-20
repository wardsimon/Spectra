function [ii,eii] = peakt(s1,bck)
%
% function [ii,eii] = peakt(s1,[nleft nright])
%
% @SPEC1D/PEAKT Trapezoidal numerical integration of peak
% in spectra s1. Returns the integrated intensity ii and
% its error ee.
%
% if specified, nleft and nright denotes how many bck points
% to take on each side.
%
% Version 2.0, February 2001
% Des McMorrow and Henrik Ronnow

%s1=struct(s1);

for il=1:length(s1)
   x=s1(il).x;
   y=s1(il).y;
   e=s1(il).e;

   y = shiftdim(y);
   e = shiftdim(e);
   m = size(y,1);

   x = x(:);
   if length(x) ~= m
       error('length(x) must equal length of first non-singleton dim of y.');
   end

   if nargin>1
      yb=sum(y([1:bck(1) end-bck(2)+1:end]))/sum(bck);
      eb=sqrt(sum(e([1:bck(1) end-bck(2)+1:end]).^2))/sum(bck);
      y=y-yb;
      e=sqrt(e.^2+eb.^2);
   end
   %   Trapezoid sum computed with vector-matrix multiply.
   ii(il)  = diff(x,1)' * (y(1:m-1) + y(2:m))/2;
   %eii= sqrt( (diff(x,1)'*e(1:m-1))^2 + (diff(x,1)'*e(2:m))^2)/2;
   eii(il) = sqrt(sum((diff(x,1).*e(1:m-1)).^2) + sum((diff(x,1).*e(2:m)).^2))/2; 
end
ii=ii(:);
eii=eii(:);

if nargout<2
   ii=[ii eii];
end

return