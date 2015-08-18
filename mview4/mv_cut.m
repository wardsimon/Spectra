function mv_cut
%
% MATLAB function to cut data from a buffer 
%
% DFM 21.4.97
%
% Modified to take monitor as a variable as well
% ARW 23.11.00

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons and validate data

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs ~=1

   mv_msg(' Must select a buffer to remove data from.');
   return

end

oper_buffs = zeros([noper_buffs,1]);
for i=1:noper_buffs

   oper_buffs(i)=str2num(get(tmv_radio(i),'String'));
   if oper_buffs(i) > nbuffs; 

      mv_msg(' The buffer does not contain data.');
      return;
 
   end

end

%----- Perform operation on buffers: cut

x  = buffers(oper_buffs(1)).xobs;
y  = buffers(oper_buffs(1)).yobs;
err= buffers(oper_buffs(1)).err;
mon= buffers(oper_buffs(1)).mon;
x_label=buffers(oper_buffs(1)).x_label;
y_label=buffers(oper_buffs(1)).y_label;
g_label=buffers(oper_buffs(1)).g_label;

mv_graph(x,y, err, x_label, y_label,g_label)

mv_msg('Click left and right limits - right button to end');

%----- Find and cut data outside of defined region

but=1;

while but==1
   
   [x_cut,y_cut,but]=ginput(2);
   if but~=1
      break
   end
   i_lt=find(x < x_cut(1));
   j_gt=find(x > x_cut(2));
   
   x_cut=[]; y_cut=[]; e_cut=[]; m_cut=[];
   if isempty(i_lt) & isempty(j_gt);
      mv_msg('    No action taken: you have removed the whole graph!')
   elseif isempty(i_lt) ~ isempty(j_gt);
      x_cut=x(j_gt(1):j_gt(length(j_gt)));
      y_cut=y(j_gt(1):j_gt(length(j_gt)));
      e_cut=err(j_gt(1):j_gt(length(j_gt)));
      m_cut=mon(j_gt(1):j_gt(length(j_gt)));	  
   elseif isempty(j_gt) ~ isempty(i_lt);
      x_cut=x(i_lt(1):i_lt(length(i_lt)));
      y_cut=y(i_lt(1):i_lt(length(i_lt)));
      e_cut=err(i_lt(1):i_lt(length(i_lt)));
      m_cut=mon(i_lt(1):i_lt(length(i_lt)));
   else
      x_cut=[x(i_lt(1):i_lt(length(i_lt))) ; ... 
          x(j_gt(1):j_gt(length(j_gt))) ];
      y_cut=[y(i_lt(1):i_lt(length(i_lt))) ; ... 
          y(j_gt(1):j_gt(length(j_gt))) ];
      e_cut=[err(i_lt(1):i_lt(length(i_lt))) ; ... 
          err(j_gt(1):j_gt(length(j_gt))) ];
      m_cut=[mon(i_lt(1):i_lt(length(i_lt))) ; ... 
          mon(j_gt(1):j_gt(length(j_gt))) ];
   end
   
   x=x_cut; y=y_cut; err=e_cut; mon=m_cut; 
   
end

g_label=['Points cut from: ' strtok(buffers(oper_buffs(1)).g_label,'.')];

%----- Update buffers

label='CUT: ';
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

%----- Plot resulta

mv_graph(x,y, err, x_label, y_label, g_label)
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return

