function mv_plot
%
% MATLAB function to plot a buffer
%
% DFM 21.4.97
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[hmv_ctrl, hmv_data]=mv_figs;
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
if length(tmv_radio)~=1

   mv_msg('Plot buffers overview');

end 
oper_buffs=str2num(char(get(tmv_radio,'String')));

m = ceil(sqrt(length(oper_buffs)));
n = ceil(length(oper_buffs)/m);

for i = 1:length(oper_buffs)

	if oper_buffs(i) > nbuffs; 

	   mv_msg(' Selected buffer does not contain data.');
	   continue; 

	end

	x  = buffers(oper_buffs(i)).xobs;
	y  = buffers(oper_buffs(i)).yobs;
	err  = buffers(oper_buffs(i)).err;
	x_label=buffers(oper_buffs(i)).x_label;
	y_label=buffers(oper_buffs(i)).y_label;
	g_label=buffers(oper_buffs(i)).g_label;

	set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',oper_buffs(i));

	%----- Plot results
	if ~isempty(hmv_data) & hmv_data
		figure(hmv_data); 
    xl=xlim;
    yl=ylim;
    subplot(m,n,1);
    xlim(xl); ylim(yl);
	end
	if i == 1, 
		mv_graph(x,y, err, x_label, y_label, g_label); 
		delete(legend);
	end
	if length(oper_buffs) > 1
		subplot(m,n,i);
		errorbar(x,y,err);
		xlabel(x_label);
		ylabel(y_label);
		title(g_label);
	end
end

%----- Reset radio buttons

mv_rtidy(0)

return
