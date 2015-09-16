function mv_combi
%
% MATLAB function to combine two buffers. If the
% x values are close than a set tolerance, then the
% points are averaged, otherwise they are appended.
%
% DFM 21.4.97
% RC 4.6.96
%
% Notes: 1. Function does not depend on order in which
%           radio buttons are pressed.
%
% Modified to take monitor as a variable as well
% ARW 22.11.00


%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons and validate data

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs < 2

   mv_msg(' Must select two or more buffers to combine.');
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

%----- Perform operation on buffers: combine

x=[]; y=[]; err=[]; mon=[];
g_label=['(' strtok(buffers(oper_buffs(1)).g_label,'.') ')'];
for i=1:noper_buffs

   buff_index=oper_buffs(i);
   x=[x ; buffers(buff_index).xobs];
   y=[y ; buffers(buff_index).yobs .* buffers(buff_index).mon];
   err=[err ; buffers(buff_index).err .* buffers(buff_index).mon];
   mon=[mon ; buffers(buff_index).mon];
   if i~=1
     g_label=['(' strtok(buffers(buff_index).g_label,'.') ')' ' + ' g_label];
   end

end

%----- Sort result in ascending order of x values

[x,perm]=sort(x);
y=y(perm);
err=err(perm);
mon=mon(perm);
min_mon=min(mon);

%----- Compare nearest neighbour points and combine if they are within the given tolerance

last_buff=buffers(oper_buffs(length(oper_buffs))).xobs;
last_buff_length=length(last_buff);

if last_buff_length > 1
   toll=abs(last_buff(last_buff_length) - ...
            last_buff(last_buff_length-1))/2;
else
   toll=abs(last_buff(last_buff_length))/10;
end

def={num2str(toll)};

tollin=inputdlg({'Input the tolerance (Xaxis):'},'Combine Tolerance',1,def);

toll=str2num(tollin{1});

xres=[];  % will contain final combined values
yres=[];
eres=[];

xcombi=[];   % will store data to be combined
ycombi=[];
ecombi=[];
mcombi=[];

for i=1:length(x);

   if isempty(xcombi)
      xcombi=[x(i)];
      ycombi=[y(i)];
      ecombi=[err(i)];
    mcombi=[mon(i)];
   else
      if x(i)-xcombi(1)>toll,
     xres=[xres;sum(xcombi)/length(xcombi)];
         yres=[yres;sum(ycombi)/sum(mcombi)];
         eres=[eres;sqrt(sum(ecombi.*ecombi))/sum(mcombi)];
         xcombi=[x(i)];
         ycombi=[y(i)];
         ecombi=[err(i)];
         mcombi=[mon(i)];
      else
         xcombi=[xcombi;x(i)];
         ycombi=[ycombi;y(i)];
         ecombi=[ecombi;err(i)];
         mcombi=[mcombi;mon(i)];
      end
   end
end;

if ~isempty(xcombi)
   xres=[xres;sum(xcombi)/length(xcombi)];
   yres=[yres;sum(ycombi)/sum(mcombi)];
   eres=[eres;sqrt(sum(ecombi.*ecombi))/sum(mcombi)];
end;

%--- Go back to variables x, y, err

x=xres;
y=yres;
err=eres;
mon=ones(size(x));
x_label=buffers(oper_buffs(1)).x_label;
y_label='Counts per monitor';

label='COMBINE: ';
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

mv_graph(x,y, err, x_label, y_label, g_label);
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return


