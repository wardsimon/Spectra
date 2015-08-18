function mf_uplot(cmd)
%
% MFIT function mf_uplot(cmd)
% Update currrent plot
% MZ 29.11.94
%

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

%=========== Initialization ====================================

if (hmf_ctrl==0) | (hmf_data==0)      % Exit if no control or data windows
   disp('Graph window not found. Load data with MFit')
   return
end
if nargin == 0, cmd=''; end
if isempty(cmd)
  cmd = 'all';
end

fitcolor  = get(findobj('tag','mf_FitColor'),'string');
if isempty(fitcolor), fitcolor = 'cyan'; end
datacolor = get(findobj('tag','mf_DataColor'),'string');
if isempty(datacolor), datacolor = 'green'; end
datacolorsel = get(findobj('tag','mf_DataColorSelected'),'string');
if isempty(datacolorsel), datacolorsel = 'magenta'; end
ebarcolor = get(findobj('tag','mf_EbarColor'),'string');
if isempty(ebarcolor), ebarcolor = 'red'; end
datline   = get(findobj('tag','mf_DataLineStyle'),'string');
if isempty(datline), datline = 'none'; end
markersize = str2num(get(findobj('Tag','mf_MarkerSize'),'String'));
if isempty(markersize)
  markersize = 10;
end

tic;
figure(hmf_data);   % Switch to data window
hax=gca;      % get current axes
set(hax,'NextPlot','Add','Units','norm','Position',[0.175 0.15 0.8 0.7]);
%set(hmf_data,'NextPlot','Add');
if ~isempty(findstr(cmd,'all'))
  cmd = [ cmd '+fit+sel+err' ];
end

limits=[get(hax,'Xlim') get(hax,'Ylim')];

%----------- Work out what to update from passed command -------

update=[0 0 0];
nnin=nargin;
if nnin>0
   if findstr(cmd,'fit')
      update(3)=1;
   end
   if findstr(cmd,'sel')
      update(2)=1;
   end
   if findstr(cmd,'err')
      update(1)=1;
   end
else
   update=[1 1 1];
end

%------------ Get handles to graph objects ---------------------

herr = findobj('Tag','mf_ebars');
hsel = findobj('Tag','mf_selected');
hfit = findobj('Tag','mf_fitline');
hfit2 = findobj('Tag','mf_fitline2');
hnsl = findobj('Tag','mf_alldat');

%============ Extract data from figure ==========================

data=get(hmf_data,'userdata');
x=data(:,1);
y=data(:,2);
err=data(:,3);
index=data(:,4);

bgcolor =get(findobj('tag','mf_FigureBgColor'),'string');
if (isempty(bgcolor) | strcmp(bgcolor, 'gray')), bgcolor = [0.8 0.8 0.8]; end
set(hmf_data,...
  'color',bgcolor,...
  'Name','MFIT: Data',...
  'visible','on',...
  'HandleVisibility','on',...
  'WindowButtonDownFcn',...
  'mf_btprs',...
  'NextPlot','Add',...
  'Interruptible','on');

%-------- Work out which values are within limits ----------------

   i = find(x>limits(1) & x<limits(2) & y>limits(3) & y<limits(4));
   if length(i)
  t = log(length(i)/100); t = max(t,.1);
  markersize = round(markersize / t);
  markersize = max(8,min(20, markersize));
   end

%   y=y(i);
%   x=x(i);
%   err=err(i);
%   index=index(i);
   if ~isempty(index)
  xi = find(index==1);
   else
  xi = [];
   end

%============= Update data ================================

   if update(1) | update(2)
      xni = find(index==0);
      if ~isempty(hnsl), delete(hnsl); end
      if ~isempty(xni)
  hnsl=line(x(xni),y(xni),...
               'Marker','.',...
               'LineStyle','none',...
               'erasemode','background',...
               'MarkerSize',ceil(markersize/2),...
               'MarkerFaceColor',datacolor,...
               'Tag','mf_alldat',...
               'color',datacolor);
      end
   end

%============= Update error bars ================================

if update(1) & ~isempty(xi)

%   refresh(hmf_data);

%------- build up nan-separated vector for ebars -----------

   mf_msg('Plotting error bars...')
   ytop = (y(xi) + err(xi))';
   ybot = (y(xi) - err(xi))';
   tee = (limits(2)-limits(1))/100;
   xleft = (x(xi)-tee)';
   xright = (x(xi)+tee)';
   nnan=NaN*ones(size(x(xi)'));

   xb=[xleft; xright; nnan; x(xi)'; x(xi)'; nnan; xleft; xright; nnan];
   yb=[ybot; ybot; nnan; ybot; ytop; nnan; ytop; ytop; nnan];
   xb=reshape(xb,9*size(xb,2),1);
   yb=reshape(yb,9*size(yb,2),1);

%---------Plot the error bars ------------------

   vis='on';
   if herr
      vis=get(herr,'visible');
   end
   h=findobj('Tag','mf_ebarsonoff');
   if strcmp(get(h,'Checked'),'on')
     vis = 'on';
   else
     vis = 'off';
   end
   if ~isempty(herr), delete(herr); end

   herr=line(xb,yb,...
             'erasemode','background',...
             'color',ebarcolor,...
             'Tag','mf_ebars',...
             'visible',vis);

end

%============ Highlight selected data points ====================

if update(2)

      if ~isempty(index)
  xi = find(index==1);
      else
  xi = [];
      end
      if ~isempty(xi)
        mf_msg('Highlighting selected data points...')
  if ~isempty(hsel), delete(hsel); end

        hsel = line(x(xi),y(xi),...
               'Marker','.',...
               'LineStyle',datline,...
               'Color',datacolorsel,...
               'erasemode','background',...
               'MarkerSize',markersize,...
               'MarkerFaceColor',datacolorsel,...
               'color',datacolorsel,...
               'Tag','mf_selected');
      else
               disp('No data selected');
      end
end

%============ Evaluate and plot fit function ===================

if update(3)

   fitfun=get(findobj('Tag','mf_FitFuncFile'),'string');
   if ~isempty(fitfun)
      tic
      p=mf_rpars;
      if ~isempty(p)
         mf_msg('Evaluating fit function...')
         h=findobj('tag','mf_FitPoints');
         if isempty(h)
            npts=max(length(x), 100);
         else
            npts = str2num(get(h,'string'));
            if isempty(npts), npts = 100;
            else npts=max(npts,length(x));
            end
         end
         npts=max(npts,length(x));
         if isempty(npts) | npts==0,  xfit=x;
         else xfit=linspace(limits(1),limits(2),npts);
         end
         xfit = xfit(:);

         % check if 1D convolution must be carried out
         hmf_convlv = findobj('Tag', 'mf_convlv');
         if ~isempty(hmf_convlv)
            UserData = get(hmf_convlv, 'UserData');
            if UserData.output.convlv & get(UserData.handles.checkbox, 'Value')
              fitfun = 'mf_convlv_func';
              mf_convlv('check',xfit); % because the number of X points may have changed !
            end
         end

         yfit=feval(fitfun,xfit,p);
         
         if ~isempty(hmf_convlv)
            UserData = get(hmf_convlv, 'UserData');
            if UserData.output.convlv & get(UserData.handles.checkbox, 'Value')
              mf_convlv('check'); %put back default x data as X axis
            end
         end

%===============================================================
         yfit = yfit(:);
         if isempty(yfit), return, end

         if length(xfit) == length(yfit)
          if ~isempty(hfit), delete(hfit); end

          hfit=line(xfit,yfit,...
                 'erasemode','background',...
                       'LineStyle','-',...
                       'color',fitcolor,...
                 'Tag','mf_fitline');

          if strcmp(fitfun,'multifunc')
            [yfit2, name, pnames, pin, separate]=multifunc(xfit, p);
            if ~isempty(hfit2), delete(hfit2); end


            for i=1:length(separate)
              if ~isempty(separate{i})

                hfit2=line(xfit,separate{i},...
                  'erasemode','background',...
                                      'LineStyle','-.',...
                                      'color','green',...
                            'Tag','mf_fitline2');
              end
            end
          end
        end
      end
   end
end
%================================================================
if isempty(herr); herr=0; end
if isempty(hsel); hsel=0; end
if isempty(hfit); hfit=0; end

set(hax,'userdata',[herr hsel hfit]);
mf_msg(['Done (' num2str(toc) 's).']);

