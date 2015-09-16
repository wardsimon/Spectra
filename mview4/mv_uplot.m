function mv_uplot(cmd)
%
% Mview function mv_uplot(cmd)
% Update currrent plot
% MZ 29.11.94
%
% mods: Now handles the error bars properly in mview DFM 18.12.95
%
[hmv_ctrl, hmv_data]=mv_figs;

%=========== Initialization ====================================
if (hmv_ctrl==0) | (hmv_data==0)			% Exit if no control or data windows
   return
end
figure(hmv_data);					% Switch to data window
hax=get(hmv_data,'CurrentAxes');			% get current axes
set(hax,'NextPlot','Add',...
        'Units','Norm',...
        'Position',[.175 .15 .8 .7]);
set(hmv_data,'NextPlot','Add');
							% Get current limits
limits=[get(hax,'Xlim') get(hax,'Ylim')];

%----------- Work out what to update from passed command -------

if nargin>0
   if strcmp(cmd,'all')
      update=[1 1 1];
   elseif strcmp(cmd,'fit')
      update=[0 0 1];
   elseif strcmp(cmd,'sel')
      update=[0 1 0];
   end
else
   update=[1 1 1];
end

%------------ Get handles to graph objects ---------------------

herronoff=findobj('Tag','mv_ebarsonoff');
hgridonoff=findobj('Tag','mv_gridonoff');
herr=findobj('Tag','mv_ebars');
hsel=max([findobj('Tag','mv_selected') 0]);
hfit=max([findobj('Tag','mv_fitline') 0]);

%============ Extract data from figure ==========================

if (update(1) | update(2))

   data=get(hmv_data,'userdata');
   [N,fourM]=size(data);                % MOD by DFM
   x=data(:,1);
   y=data(:,2);
   std=data(:,3);
   index=data(:,4);

%-------- Work out which values are within limits 

   i= x>limits(1) & x<limits(2) & y>limits(3) & y<limits(4);
   y=y(i);
   x=x(i);
   std=std(i);
   index=index(i);

end

%============= Update error bars ================================

if update(1)

   refresh(hmv_data);

%----- Loop over the number of scans added by DFM

   xb_sav=[]; yb_sav=[];
   for isc=1:fourM/4

      x  =data(:,1+(isc-1)*4); x  =  x(~isnan(x  ));
      y  =data(:,2+(isc-1)*4); y  =  y(~isnan(y  ));
      std=data(:,3+(isc-1)*4); std=std(~isnan(std));
	
%------- build up nan-separated vector for ebars -----------

      ytop = (y + std)';
      ybot = (y - std)';
      tee = (limits(2)-limits(1))/100;
   
      xleft = (x-tee)';
      xright = (x+tee)';
      nnan=NaN*ones(size(x'));

      xb=[xleft; xright; nnan; x'; x'; nnan; xleft; xright; nnan];
      yb=[ybot; ybot; nnan; ybot; ytop; nnan; ytop; ytop; nnan];
      xb=reshape(xb,9*size(xb,2),1);
      yb=reshape(yb,9*size(yb,2),1);
      xb_sav=[xb_sav ; xb];
      yb_sav=[yb_sav ; yb];

   end

%---------Plot the error bars ------------------

%------------ Check whether grid and errorbars are desired ------
   viserr='on';
   if herr
      viserr=get(herr,'visible');
      delete(herr);
   elseif strcmp(get(herronoff,'checked'),'off')
      viserr='off';
   end
   figure(hmv_data);	
   ebarcolor=get(findobj('Tag','mv_EbarColor'),'String');
   herr=line(xb_sav,yb_sav,'erasemode','background',...
                           'color',ebarcolor,'Tag','mv_ebars',...
                           'visible',viserr);
                        
   visgrid='on';
   if strcmp(get(hgridonoff,'checked'),'off')
      visgrid='off';
   end;
   set(gca,'Xgrid',visgrid);
   set(gca,'Ygrid',visgrid);
                        
end


%============ Highlight selected data points ====================

if update(2)

   x=data(:,1);
   y=data(:,2);
   std=data(:,3);
   index=data(:,4);

%-------- Work out which values are within limits ----------------

   i= x>limits(1) & x<limits(2) & y>limits(3) & y<limits(4);
   y=y(i);
   x=x(i);
   std=std(i);
   index=index(i);

   if hsel delete(hsel); end
   figure(hmv_data);
   datacolor=get(findobj('Tag','mv_DataColor'),'String');	
   Marker_Size=str2num(get(findobj('Tag','mv_MarkerSize'),'String'));
   hsel=line(x(find(index==1)),y(find(index==1)),...
                               'Marker','o',...
                               'LineStyle','none',...
                               'erasemode','background',...
                               'MarkerSize',Marker_Size,...
                               'MarkerFaceColor',datacolor,...
                               'color',datacolor,...
                               'Tag','mv_selected');
end

%============ Evaluate and plot fit function ===================

if update(3)
%h=get(findobj('Tag','hmv_file'),'Userdata');
   fitfun='';
   fundir='';

   if ~isempty(fitfun)
      p=mv_rpars;
      if p~=[]
         xfit=limits(1):(limits(2)-limits(1))/100:limits(2);
         yfit=feval(fitfun,xfit,p);
            if hfit delete(hfit); end
               figure(hmv_data);	
               hfit=line(xfit,yfit,'erasemode','background','color','c',...
	    			   'Tag','mv_fitline');
	    end
    end
end

if isempty(herr); herr=0; end
if isempty(hsel); hsel=0; end
if isempty(hfit); hfit=0; end

set(hax,'userdata',[herr hsel hfit]);
refresh


