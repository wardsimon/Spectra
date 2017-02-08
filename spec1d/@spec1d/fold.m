function svec=fold(svec, smodell)
% function s=smooth(s, smodell)
%
%fragment
%
% SPED1D/Smoothes spectrum s with smodell
% Performs a convolution with smodell
%
% CK 20.01.08

%----- Get bin centers from dx
mx=smodell.x;
my=smodell.y;

for nspec=1:length(svec)
   
   s=svec(nspec);
   nps=length(s.x);
   y=[]; e=[];
   rangex=max(s.x)-min(s.x);
   if rangex>max(mx)
        mx=[mx;rangex];
        my=[my;0];
        [mx,ind]=sort(mx);
        my=my(ind);
   elseif -rangex<min(mx)
        mx=[mx;-rangex];
        my=[my;0];
        [mx,ind]=sort(mx);
        my=my(ind);
   end

   for n=1:nps
      %r=exp(-(s.x-s.x(n)).^2/dx^2/2);
      r=interp1(mx,my,s.x-s.x(n));
      r=r/sum(r); 
      y(n)=r'*s.y;
      e(n)=sqrt(sum((s.e.*r).^2));
   end

   s.y=y(:);
   s.e=e(:);

   svec(nspec)=s;
   
end

return