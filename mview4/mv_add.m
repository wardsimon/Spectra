function mv_add
%
% MATLAB function to add two or more buffers.
% Also accounts for monitor
% Adapted from mv_combi
% ARW 22.11.00
%
% mv_combi written by
% DFM 21.4.97
% RC 4.6.96
%
% Notes: 1. Function does not depend on order in which
%           radio buttons are pressed.
%


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

%----- Perform operation on buffers: add

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

%----- Check to see if scans overlap
x_lims=[min(buffers(oper_buffs(1)).xobs) max(buffers(oper_buffs(1)).xobs)];
for i=2:noper_buffs
  if x_lims(1) > max(buffers(oper_buffs(i)).xobs) | ...
     x_lims(2) < min(buffers(oper_buffs(i)).xobs)
    mv_msg(' Scans don''t overlap');
    return
  end
  if min(buffers(oper_buffs(i)).xobs) > x_lims(1)
    x_lims(1) = min(buffers(oper_buffs(i)).xobs);
  end;
  if max(buffers(oper_buffs(i)).xobs) < x_lims(2)
    x_lims(2) = max(buffers(oper_buffs(i)).xobs);
  end
end

%----- Choose data within limits and chuck the rest
i = find(x>=x_lims(1) & x<=x_lims(2));
x=x(i); y=y(i); err=err(i); mon=mon(i);


%----- Compare nearest neighbour points and add if they are within the given tolerance

last_buff=buffers(oper_buffs(length(oper_buffs))).xobs;
last_buff_length=length(last_buff);

if last_buff_length > 1
   toll=abs(last_buff(last_buff_length) - ...
            last_buff(last_buff_length-1))/2;
else
   toll=abs(last_buff(last_buff_length))/10;
end

def={num2str(toll)};

tollin=inputdlg({'Input the tolerance (X axis):'},'Add Tolerance',1,def);

toll=str2num(tollin{1});

xres=[];  % will contain final added values
yres=[];
eres=[];
mres=[];

xcombi=[];   % will store data to be added
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
         yres=[yres;sum(ycombi)];
         eres=[eres;sqrt(sum(ecombi.*ecombi))];
     mres=[mres;sum(mcombi)];
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
   yres=[yres;sum(ycombi)];
   eres=[eres;sqrt(sum(ecombi.*ecombi))];
   mres=[mres;sum(mcombi)];
end;

%--- Go back to variables x, y, err, mon

x=xres;
y=yres;
err=eres;
if mean(mres) == noper_buffs
  mon=ones(size(y));
  y_label='Counts';
else
  mon=mres;
  y_label='Counts per monitor';
end;
x_label=buffers(oper_buffs(1)).x_label;

y = y./mon;
err = err./mon;

%----- Update buffers

label='ADD: ';
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

mv_graph(x,y, err, x_label, y_label, g_label);
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return
