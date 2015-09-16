function [partext,h]=formatpars(fpars,flag,press)
%
% function partext=formatpars(fpars,flag)
%
% @SPEC1D/FORMATPAR Writes parameter names and values to string partext  
%
%                   If flag>0, write to screen.
%                      flag=1, Top right corner.
%                      flag=2, Top left corner.
%                      flag=3, Lower left corner.
%                      flag=4, Lower right corner.
%
% DFM 1.10.98
%

if nargin<2
   flag=1;
end

if nargin<3
   press=2;
end

[npars,dummy]=size(fpars.pvals);
s1=[];
for i=1:npars, s1=strvcat(s1,' = '); end
s2=[];
%for i=1:npars, s2=strvcat(s2,' +/- '); end
for i=1:npars, s2=strvcat(s2,' \pm '); end

%partext=[strvcat(fpars.pnames) s1 strvcat(num2str(fpars.pvals,['%0.',num2str(press+1),'g'])) s2 ...
%                           strvcat(num2str(fpars.evals,['%0.',num2str(press-1),'g']))];
partext=[strvcat(fpars.pnames) s1 strvcat(num2str(fpars.pvals,['%0.',num2str(press),'f'])) s2 ...
                           strvcat(num2str(fpars.evals,['%0.',num2str(press),'f']))];

if flag>0

 s1d_h=findobj('Tag','s1d_DataWindow');
 if isempty(s1d_h), 
%      disp(' No data window, using current!')
    s1d_h=gcf;
 end
 s1d_child=get(s1d_h,'Children');
 axes(s1d_child(1))
 if flag==1
   h=text(0.95,0.95,partext,'units','normalized','horizontalalignment','right','verticalalignment','top');
 elseif flag==2
   h=text(0.05,0.95,partext,'units','normalized','horizontalalignment','left','verticalalignment','top');
 elseif flag==3
   h=text(0.05,0.05,partext,'units','normalized','horizontalalignment','left','verticalalignment','bottom');
 elseif flag==4
   h=text(0.95,0.05,partext,'units','normalized','horizontalalignment','right','verticalalignment','bottom');
 else
   h=text(0.95,0.95,partext,'units','normalized','horizontalalignment','right','verticalalignment','top');
 end
end
