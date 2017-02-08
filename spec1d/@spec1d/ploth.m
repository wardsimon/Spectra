function [hout,hbout]=plot(varargin)
%
% function plot(varargin)
%
% SPEC1D/PLOT Plot for spectra
%
% Usage: 1. Simple plot or overlay of spectra: >> plot(s1,s2,s3)
%                                              >> plot(s1,10.*s2)
%        2. Possibility to choose axes type:   >> plot(s1,'semilogy')
%                                              >> plot(s1,'semilogx')
%                                              >> plot(s1,'loglog')
%
% Tags:  Data window: s1d_DataWindow 
%
% DFM 2.10.98

% HMR 21.5.99 Outputs handle to markers and to errorbars (if asked for).
r.x=[];
r.y=[];
r.e=[];

%----- Draw data window

%plot_dwin;

oldvarargin=varargin;
varargin={};
for i=1:length(oldvarargin)
   if isa(oldvarargin{i},'spec1d')
     for n=1:length(oldvarargin{i})
       varargin{end+1}=oldvarargin{i}(n);
     end
   else
       varargin{end+1}=oldvarargin{i};   
   end  
end

%----- Work out limits on graph
for i=1:length(varargin)
   if isa(varargin{i},'spec1d')
     r.x=[r.x; varargin{i}.x];
     r.y=[r.y; varargin{i}.y];
     r.e=[r.e; varargin{i}.e];
   end  
end
[r.x,sort]=sort(r.x);
r.y=r.y(sort);
r.e=r.e(sort);
maxx=max(r.x); minx=min(r.x); xr=maxx-minx;
maxy=max(r.y); miny=min(r.y); yr=maxy-miny;
maxe=max(r.e); mine=min(r.e); er=maxe-mine;

%axis([minx-xr*0.1 maxx+xr*0.1 miny-(yr+er)*0.1 maxy+(yr+er)*0.1]);
%axis([miny-(yr+er)*0.1 maxy+(yr+er)*0.1 minx-xr*0.1-eps maxx+xr*0.1+eps]);
hold on
colororder=[0 0 1; 0 1 0; 1 0 0; 0 0 0; 1 0 1; 0 1 1];
markerorder=str2mat('o','square','diamond','^','>','<');


hout=[];
hbout=[];
for i=1:length(varargin)

   if isa(varargin{i},'spec1d')
   
%----- Plot graph and error bars

      x=varargin{i}.x; y=varargin{i}.y; err=varargin{i}.e; yfit=varargin{i}.yfit;
      x=x(:); y=y(:); err=err(:); yfit=yfit(:);

%----- build up nan-separated vector for bars

      ytop = (y + err)';
      ybot = (y - err)';
      tee = (max(x)-min(x))/100;  % make tee .02 x-distance for error bars
%      tee = 0;
      xleft = (x-tee)';
      xright = (x+tee)';
      nnan=NaN*ones(size(x'));

      xb=[xleft; xright; nnan; x'; x'; nnan; xleft; xright; nnan];
      yb=[ybot; ybot; nnan; ybot; ytop; nnan; ytop; ytop; nnan];
      n=9*size(xb,2);

      xb=reshape(xb,n,1);
      yb=reshape(yb,n,1);
   
%----- Draw lines
   
      hle=line(yb,xb,'color',[.5 .5 .5]);
      set(hle,'ButtonDownFcn','editline(gco);');
      hll=line(y,x,'color',colororder(1+mod(i-1,6),:),...
         'LineStyle','None','Marker',deblank(markerorder(1+mod(i-1,6),:)),...
         'MarkerSize',7,'MarkerFaceColor',colororder(1+mod(i-1,6),:));
      if mod(i,2)==0
         set(hll,'markerfacecolor','none')
%         set(hll,'markerfacecolor','w')
      end
      set(hll,'ButtonDownFcn','editline(gco);');
   
      if ~isempty(yfit); 
         hlf=line(yfit,x,'color','m','LineStyle','-','LineWidth',2,'Marker','none');
         set(hlf,'ButtonDownFcn','editline(gco);');
      end
   end
   hout=[hout;hll];   
   hbout=[hbout;hle];   
end
if nargout<1
   clear hout hbout
elseif nargout<2
   clear hbout
end

set(gca,'Box','on')

%----- Make title

set(get(gca,'Title'),...              
   'String',varargin{1}.datafile,...
   'Tag','plot_text_title',...
   'FontName','Times',...
   'Fontsize',14,...
   'buttondownfcn','edtext',...
   'handlevisibility','on',...
   'interpreter','none');

set(get(gca,'Xlab'),...              
   'String',varargin{1}.x_label,...
   'Tag','plot_text_xlabel',...
   'FontName','Times',...
   'Fontsize',14,...
   'buttondownfcn','edtext',...
   'handlevisibility','on',...
   'interpreter','none');

set(get(gca,'Ylab'),...              
   'String',varargin{1}.y_label,...
   'Tag','plot_text_ylabel',...
   'FontName','Times',...
   'Fontsize',14,...
   'buttondownfcn','edtext',...
   'handlevisibility','on',...
   'interpreter','none');

set(gca,'FontName','Times','FontSize',14);
set(gca,'Xgrid','on','Ygrid','on');
set(gca,'ButtonDownFcn','scribeaxesdlg(gca);');

%----- Check for log axes

axes_choice=strvcat('semilogy','semilogx','loglog');

for i=1:length(varargin)

   if isa(varargin{i},'char')

      axes_type=varargin{i};
      switch axes_type
         case 'semilogy', set(gca,'Yscale','log');
         case 'semilogx', set(gca,'Xscale','log');
         case 'loglog'  , set(gca,'Xscale','log','Yscale','log');
         otherwise, disp('Not a valid axis type')
      end

   end      

end

%----- Start of subroutines

%======================================================================================
function plot_dwin
%
% MATLAB function to draw data window and add menus
%
%
%
hplot=findobj('Tag','s1d_DataWindow');

if ~isempty(hplot)
   figure(hplot)
   cla
else
   hplot=figure('Position',[50 50 500 400],...
          'Tag','s1d_DataWindow',...
          'Visible','on',...
          'Name','Spectra Plot',...
          'Windowbuttondownfcn','edtext(''hide'')',...         
          'Menubar','None',...
          'NumberTitle','off');

%---------- Print menu-----------------------------

   hprt=uimenu(hplot,'Label','Print');
   uimenu(hprt,'Label','Print figure',...
            'Callback',['hplot=findobj(''Tag'',''s1d_DataWindow'');'...
                        'figure(hplot);'...
                        'set(gca,''units'',''points'');'...
                        'pos=get(gca,''Position'');'...
                        'set(gcf,''papertype'',''a4letter'');'...
                        'set(gcf,''units'',''centimeters'');'...
                        'set(gcf,''paperunits'',''centimeters'');'...
                        'p=get(gcf,''position'');'...
                        'width=p(3); height=p(4);'...
                        'paper=get(gcf,''papersize'');'...
                        'ppos=[(paper(1)-width)/2 (paper(2)-height)/1.25 width height];'...
                        'set(gcf,''paperposition'',ppos);'...
                        'printdlg(gcf);']);
                       
   uimenu(hprt,'Label','Save figure',...
      'Callback','[f,p] = uiputfile(''*.m'',''Save Figure As...'');if f,eval([''print -dmfile '' p f ]); disp(''Figure saved''); end');
       
end