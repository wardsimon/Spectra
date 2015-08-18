function mv_intp
%
% MATLAB function to perform interpolation on two buffers.
% in mview.
%
%
% DFM 26.3.98
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs ~= 2

   mv_msg(' Must select two buffers only for this operation.');
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

%----- Perform one of the allowed operations
x_length_test = length(buffers(oper_buffs(1)).xobs);

y=zeros([length(x_length_test),1]);              
err=zeros([length(x_length_test),1]);
g_label=[];

x1  = buffers(oper_buffs(1)).xobs;
x2  = buffers(oper_buffs(2)).xobs;
y1  = buffers(oper_buffs(1)).yobs;
y2  = buffers(oper_buffs(2)).yobs;
e1  = buffers(oper_buffs(1)).err;
e2  = buffers(oper_buffs(2)).err;
m1  = buffers(oper_buffs(1)).mon;
m2  = buffers(oper_buffs(2)).mon;
x_label_1=buffers(oper_buffs(1)).x_label; 
x_label_2=buffers(oper_buffs(2)).x_label;
y_label_1=buffers(oper_buffs(1)).y_label;
y_label_2=buffers(oper_buffs(2)).y_label;
g_label_1=strtok(buffers(oper_buffs(1)).g_label,'.');
g_label_2=strtok(buffers(oper_buffs(2)).g_label,'.');

if length(x1)==length(x2) 
   if isempty(find(x1~=x2))
      mv_msg(' No interpolation:  X ranges are identical');
      return
   end
end

x2_length=length(x2);

while x2(1) < min(x1); x2(1)=[]; y2(1)=[]; e2(1)=[]; m2(1)=[]; end      
while x2(end) > max(x1); 
   x2(end)=[]; 
   y2(end)=[]; 
   e2(end)=[];
   m2(end)=[];
end
x     = x2;
y     = interp1(x1,y1,x2);
err   = interp1(x1,e1,x2);
mon	  = interp1(x1,m1,x2);
x_label = x_label_1;
y_label = 'Counts';
g_label = ['Interpolated: ' '(' g_label_1 ')' ' on ' '(' g_label_2 ')'];
label='INTER: ';

%----- First return x2 buffer if it has been truncated
 
if length(x2) ~= x2_length
       
   [incr_buffs]=mv_ubuff(x2,y2,e2,m2,x_label_2,y_label_2,g_label_2,[]);
   set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',incr_buffs);

%----- Change text in window to file name and plot results

   hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
   label2=['Truncated buffer' num2str(oper_buffs(2)) '  returned'];
   set(hmv_text(incr_buffs),'String',label2);
   
end   

[incr_buffs]=mv_ubuff(x,y,err,mon,x_label,y_label,g_label,[]);
set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',incr_buffs);

%----- Change text in window to file name and plot results

hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
for i=1:noper_buffs; label=[label, num2str(oper_buffs(i)) '  ']; end
set(hmv_text(incr_buffs),'String',label,'TooltipString',label);

%----- Plot results

mv_graph(x,y,err,x_label,y_label,g_label);
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return	
