function stats=peakm(s1,bg_pts)
%
% function stats=peakm(s1,[nleft nright])
%
% SPEC1D/PEAK Calculates properties of peak in s1
%             from its moments. 
%             s1 may be an array of spec1d objects.
%
% Notes:
%
% bg_pts=[nleft nright] and is the number of points in
% the left and right hand wings of the peak
% used in calculating the background. If not
% specified no background correcgtion is made.
%
% stats(1)= I_cor  : Peak intensity corrected for background
% stats(2)= Centre : First moment of peak
% stats(3)= FWHM   : Second moment of peak
% stats(4)= I_int  : Integrated area of peak= I_cor*FWHM
% stats(5)= Bg     : Background
% stats(6)= I_aver : Average Y
% stats(7)= I_std  : Standard deviation in Y
%
% Version 2.0, February 2001
% Des McMorrow and Henrik Ronnow

%s1=struct(s1);
for il=1:length(s1)

   x=s1(il).x; y=s1(il).y; e=s1(il).e;
   npts=length(x);

%----- Calculate average background

   if nargin == 1
   
      Bg=0;
   
   elseif nargin==2 & length(bg_pts) <= 2
   
      nleft=bg_pts(1);
      nright=bg_pts(end);
   
      Bg_left=0; Bg_right=0;
   
      if nleft~=0;  Bg_left =sum(y(1:nleft))/nleft; end
      if nright~=0; Bg_right=sum(y(end-nright+1:end))/nright; end
   
      if nleft~=0 & nright ~=0
         Bg=(Bg_left+Bg_right)/2; 
      elseif nleft~=0 & nright ==0
         Bg=Bg_left;
      elseif nleft==0 & nright ~=0
         Bg=Bg_right;
      end
   end   

   y=y-Bg;

%----- Store results in output array

   stats(il,1)=max(y);
   stats(il,2)=sum(x.*y)/sum(y);
   stats(il,3)=sqrt(abs(sum(x.^2.*y)/sum(y)-(sum(x.*y)/sum(y))^2));
   stats(il,4)=stats(il,1)*stats(il,3); 
   stats(il,5)=Bg;
   stats(il,6)=mean(y);
   stats(il,7)=std(y);

end
