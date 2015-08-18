function sout=interpolate(s, xnew)
% function s=interpolate(s, xnew)
%
% SPED1D/Interpolates a spectrum to new x-values
%
% At the moment only does linear interpolation.
% Could be extended to quadratic and other methods.
%
% HMR, NBC 20.11.2000
%

if nargin<2
  error('Not enough input')
end

if isa(xnew,'spec1d')
  xnew=getfield(xnew,'x');
elseif length(xnew)==1
  xnew=min(s.x):xnew:max(s.x);
elseif isa(xnew,'char') & strcmp(xnew,'auto')
  xnew=2*mean(s.x(2:end)-s.x(1:end-1));  
  xnew=min(s.x):xnew:max(s.x);
end

x=s.x;
y=s.y;
e=s.e;

for j=1:length(xnew)
   n1=find(x<=xnew(j));
   n2=find(x> xnew(j));
   if isempty(n1)
      n1=n2(1);
      n2=n2(2);
   elseif isempty(n2)
      n2=n1(end);
      n1=n1(end-1);
   else
      n1=n1(end);
      n2=n2(1);
   end
   s.y(j)=y(n1)*(x(n2)-xnew(j))/(x(n2)-x(n1))+...
          y(n2)*(x(n1)-xnew(j))/(x(n1)-x(n2));   
   s.e(j)=sqrt((e(n1)*(x(n2)-xnew(j))/(x(n2)-x(n1)))^2+...
               (e(n2)*(x(n1)-xnew(j))/(x(n1)-x(n2)))^2);
end

sout = spec1d(xnew,s.y(1:length(xnew)),s.e(1:length(xnew)));