function [incr_buffs]=mv_ubuff(x,y,err,mon,x_label,y_label,g_label,datafile)
%
% MVIEW function to update the data stored in buffers
%
% DFM 21.4.97
%
% Modified to take monitor as a variable as well
% ARW 22.11.00
%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');

%----- Update buffers

buffers_local.xobs=x;
buffers_local.yobs=y;
buffers_local.err=err;
buffers_local.mon=mon;
buffers_local.x_label=x_label;
buffers_local.y_label=y_label;
buffers_local.g_label=g_label;
if isempty(datafile), datafile = g_label; end
buffers_local.datafile=datafile;

if isempty(buffers)
   buffers=buffers_local;
else
   buffers=[buffers buffers_local];
end

[dummy,incr_buffs]=size(buffers);

set(findobj('Tag','hmv_buffers'),'Userdata',buffers);
