function mv_print
%
% MATLAB function to print a data window in a tidy way
%
%
% DFM 30.5.97
%
hmv_DataWindow=findobj('Tag','MVIEW: Data');
figure(hmv_DataWindow)
set(gca,'units','points')
pos=get(gca,'Position');

hmv_AllData=findobj('Tag','mv_AllData');
hmv_AllDatastr=findobj('Tag','mv_AllDatastr');

fontsize=12;
if ~isempty(hmv_AllData) | ~isempty(hmv_AllDatastr)

   Data=get(hmv_AllData,'Userdata');
   Datastr=get(hmv_AllDatastr,'Userdata');

   minData=min(Data);
   maxData=max(Data);
   meanData=mean(Data);
   stdData=std(Data);

   [npts,nvars]=size(Data);
   
%   s=sprintf(' \t   Min \t \t   Max \t \t   Mean \t  Std. Dev. ');
   for i=1:nvars
 
%      s=sprintf(' %s :\t %0.3e \t %0.3e \t %0.3e \t %0.3e  ', ...
%                Datastr(i,:),minData(i),maxData(i),meanData(i),stdData(i));

%      text(0,0,s,'units','points','position',...
%           [pos(1)-50 pos(2)+100-(i-1)*fontsize*1.5],'fontsize',fontsize)

   end

end

set(gcf,'papertype','a4letter');
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
p=get(gcf,'position');
width=p(3); height=p(4);

paper=get(gcf,'papersize');
ppos=[(paper(1)-width)/2 (paper(2)-height)/1.25 width height];
set(gcf,'paperposition',ppos);

%if strcmp(computer,'PCWIN')
%   print -dwin -v
%else
%   print
%end

printdlg(gcf);

