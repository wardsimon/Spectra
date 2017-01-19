function svec=smooth(svec, dx)
% function s=smooth(s, dx)
%
% SPED1D/Smoothes spectrum s with gaussian width dx
% Performs a convolution with gaussian of width dx
%
% HMR 01.10.99

%----- Get bin centers from dx

for nspec=1:length(svec)
   
   s=svec(nspec);

   nps=length(s.x);

   if nargin<2 % guess a reasonable smooting width.
      dx=mean(diff(s.x));
   end

   if length(dx)~=1
     error('Smoothing width must be scalar')   
   end

   y=[]; e=[];
   for n=1:nps
      r=exp(-(s.x-s.x(n)).^2/dx^2/2);
      r=r/sum(r); 
      y(n)=r'*s.y;
      e(n)=sqrt(sum((s.e.*r).^2));
   end

   s.y=y(:);
   s.e=e(:);

   svec(nspec)=s;
   
end

return