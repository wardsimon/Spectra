function svec=cut(svec,varargin)
%
% function r=cut(s1,[xleft1 xright1],[xleft2 xright2],...)
%
% @SPEC1D/CUT Cuts data from a spec1d spectrum or spec1d array 
%             using the xranges specified by [xleft1 xright1], etc.
%
% [NaN   xright] : cuts data greater than xright
% [xleft NaN   ] : cuts data less than xleft1
% [xleft xright] : If xleft < xright, cuts data < xleft and > xright
% [xleft xright] : If xleft > xright, cuts data > xleft and < xright
%
% Example:
%
% Cut data outside of the range 101-102
% >r=cut(s1,[101 102]);
% Cut data in the x range 101-102 and 109-110 from s1
% >r=cut(s1,[102 101],[110 109]);
%
% Version 2.0, April 2001
% Des McMorrow and Henrik R�nnow
% speed increased by HMR 2003/3

for n=1:length(svec)
   
   if nargin ==1
      disp('Spec1d Error: Must specify range to cut')
      return
   end
   if length(varargin{1}) ~= 2
      disp('Spec1d Error: Must specify cut range by two variables')
      return
   end

   for i=1:length(varargin)

      xleft=varargin{i}(1);
      xright=varargin{i}(2);
   
      if xleft > xright
         xci=find(svec(n).x<xleft&svec(n).x>xright);
         svec(n).x(xci)=[]; 
         svec(n).y(xci)=[]; 
         svec(n).e(xci)=[]; 
         if ~isempty(svec(n).yfit)
             svec(n).yfit(xci)=[]; 
         end
      else
         xcl=find(svec(n).x<xleft);
         xcr=find(svec(n).x>xright);
         svec(n).x([xcl;xcr])=[]; 
         svec(n).y([xcl;xcr])=[]; 
         svec(n).e([xcl;xcr])=[]; 
         if ~isempty(svec(n).yfit)
             svec(n).yfit([xcl;xcr])=[]; 
         end
      end
   
   end   

end