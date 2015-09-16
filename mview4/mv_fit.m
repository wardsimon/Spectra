function mv_fit
%
% MATLAB function to fit buffer using mfit
%
% DFM 21.4.97
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs~=1

   mv_msg(' Must select only one buffer to Fit.');
   return

end 

oper_buffs = zeros([noper_buffs,1]);
for i=1:noper_buffs

   oper_buffs(i)=str2num(get(tmv_radio(i),'String'));
   if oper_buffs(i) > nbuffs; 

      mv_msg(' Selected buffer does not contain data.');
      return; 

   end

end

%----- Perform operation on buffers: fit

buff_fit = oper_buffs;

x=buffers(buff_fit).xobs;
y=buffers(buff_fit).yobs;
err=buffers(buff_fit).err;
xlab=buffers(buff_fit).x_label;
ylab=buffers(buff_fit).y_label;
name = buffers(buff_fit).g_label;

mfitgo(x,y,err,xlab,ylab,name);
   
%----- Reset radio buttons

mv_rtidy(0);

return

