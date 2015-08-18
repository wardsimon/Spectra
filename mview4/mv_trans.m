function mv_trans(opera)
%
% MATLAB function to transform a buffer.
%
% DFM 21.4.97
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs ~=1

   mv_msg(' Must select only one buffer to transform.');
   return

end

oper_buffs = zeros([noper_buffs,1]);
for i=1:noper_buffs

   oper_buffs(i)=str2num(get(tmv_radio(i),'String'));
   if oper_buffs(i) > nbuffs; 

      mv_msg(' Not all buffers contain data.');
      return; 
   
   end

end

%----- Perform the operation: transform

x  = buffers(oper_buffs(1)).xobs;
y  = buffers(oper_buffs(1)).yobs;
err= buffers(oper_buffs(1)).err;
mon= buffers(oper_buffs(1)).mon;
x_label=buffers(oper_buffs(1)).x_label;
y_label=buffers(oper_buffs(1)).y_label;
g_label=buffers(oper_buffs(1)).g_label;

eval([opera ';']);
label='TRANS: ';
tooltip=label;
for i=1:noper_buffs; 
	label=[label, num2str(oper_buffs(i)) '  ']; 
	tooltip = [tooltip sprintf('\n') buffers(oper_buffs(i)).datafile];
end

%----- Update buffers

[incr_buffs]=mv_ubuff(x,y,err,mon,x_label,y_label,g_label,[]);
set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',incr_buffs);

%----- Change text in window to file name and plot results

hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
set(hmv_text(incr_buffs),'String',label,'TooltipString',tooltip);

%----- Plot results

mv_graph(x,y, err, x_label, y_label, g_label)
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return

