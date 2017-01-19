function rvec=cuty(svec,varargin)
%
% function r=cuty(s1,[ybot1 ytop1],[ybot2 ytop2],...)
%
% @SPEC1D/CUT Cuts data from a spec1d spectrum or spec1d array 
%             using the yranges specified by [ybot1 xtop1], etc.
%
% [NaN   ytop]   : cuts data greater than ytop
% [ybot  NaN ]   : cuts data less than xbot
% [ybot ytop] : If ybot < ytop, cuts data < ybot and > ytop
% [ybot ytop] : If ybot > ytop, cuts data > ybot and < ytop
%
% Example:
%
% Version 2.0, May 2002
% Des McMorrow and Henrik Rønnow

for ns=1:length(svec)
   
   s1=svec(ns);
   
   x=s1.x; y=s1.y; e=s1.e; yfit=s1.yfit;
   x_label=s1.x_label; y_label=s1.y_label;

   if nargin ==1
      disp('Spec1d Error: Must specify range to cut')
      return
   end
   if length(varargin{1}) ~= 2
      disp('Spec1d Error: Must specify cut range by two variables')
      return
   end

   for i=1:length(varargin)

      ybot=varargin{i}(1);
      ytop=varargin{i}(2);
   
      if ybot > ytop
         yci=find(y<ybot&y>ytop);
         x(yci)=[]; y(yci)=[]; e(yci)=[]; 
         if ~isempty(yfit), yfit(yci)=[]; end
      else
         ycb=find(y<ybot);
         yct=find(y>ytop);
         x([ycb;yct])=[]; y([ycb;yct])=[]; e([ycb;yct])=[]; 
         if ~isempty(yfit), yfit([ycb;yct])=[]; end
      end
   
   end   

   r.x=x;
   r.y=y;
   r.e=e;
   r.x_label=x_label;
   r.y_label=y_label;
   r.datafile=s1.datafile;
   r.yfit=yfit;
   r=spec1d(r);
   rvec(ns)=r;
   
end