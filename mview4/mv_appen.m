function mv_appen
%
% MATLAB function to append two or more buffers together
% in mview.
%
% DFM 21.4.97
%
% Notes: 1. Function does not depend on order in which
%           radio buttons are pressed.

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs < 2

   mv_msg(' Must select two or more buffers to append.');
   return

end 

oper_buffs = zeros([noper_buffs,1]);
for i=1:noper_buffs

   oper_buffs(i)=str2num(get(tmv_radio(i),'String'));
   if oper_buffs(i) > nbuffs; 

      mv_msg(' Not all selected buffers contain data.');
      return; 
   
   end

end

%----- Sort radio buttons

oper_buffs=sort(oper_buffs);

%----- Perform operation on buffers: append

x=[]; y=[] ; err=[]; mon=[]; 

g_label=['(' strtok(buffers(oper_buffs(1)).g_label,'.') ')'];
for i=1:noper_buffs

    buff_index=oper_buffs(i);
    x  = [x  ; buffers(buff_index).xobs];
    y  = [y  ; buffers(buff_index).yobs];
    err= [err; buffers(buff_index).err];
    mon= [mon; buffers(buff_index).mon];
    if i~=1
       g_label=['(' strtok(buffers(buff_index).g_label,'.') ')' ' + ' g_label];
    end

end
    
x_label=buffers(oper_buffs(1)).x_label;
y_label='Counts';

%----- Sort data according to the value of x

[x,idx]=sort(x);
y=y(idx);
err=err(idx);
mon=mon(idx);

%----- Update buffers

label='APPEND: ';
tooltip = label;
tooltip = label;
for i=1:noper_buffs; 
	label=[label, num2str(oper_buffs(i)) '  ']; 
	tooltip = [tooltip sprintf('\n') buffers(oper_buffs(i)).datafile];
end

[incr_buffs]=mv_ubuff(x,y,err,mon,x_label,y_label,g_label,[]);
set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',incr_buffs);

%----- Change text in window to file name and plot results

hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
set(hmv_text(incr_buffs),'String',label,'TooltipString',tooltip);

%----- Plot results

mv_graph(x,y, err,x_label, y_label,g_label);
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return
