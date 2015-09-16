function mv_graph(x,y, err, x_label, y_label, g_label)
%
% MVIEW function to plot data
%
% DFM 12.4.95
%
hmv_data=mv_dwin(x_label, y_label, g_label);
figure(hmv_data);

%---------- Attach data to userdata ------------------

set(hmv_data,'userdata',[x y err ones(size(y))]);

%----------------- Set limits  --------------------

h=get(findobj('Tag','mv_fixax'),'Checked');

if strcmp(h,'off')
   yrange=max(y)-min(y);
   if yrange==0
      yrange=yrange+max([1e-3*mean(y) 1e-6]);
   end
   xrange=max(x)-min(x);
   if xrange==0
      xrange=xrange+max([1e-3*mean(x) 1e-6]);
   end
   x_lim=[min(x)-0.01*xrange max(x)+0.01*xrange];
   y_lim=[min(y)-0.1*yrange max(y)+0.1*yrange];
else
   x_lim=get(gca,'Xlim');
   y_lim=get(gca,'Ylim');
end

set(gca,'Xlim',x_lim);
set(gca,'Ylim',y_lim);

%----------- Do the plot ------------------------------

mv_uplot;
