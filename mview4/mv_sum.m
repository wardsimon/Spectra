function mv_sum

%  Function to add the data and produce an integrated intensity.
% Does this by the trapezoid rule, where the data points are
% joined by a line, and the areas of the trapezoids are calculated.
% 
% ARW 12-10-97

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
x_label=buffers(oper_buffs(1)).x_label;
y_label=buffers(oper_buffs(1)).y_label;
g_label=buffers(oper_buffs(1)).g_label;

x_s = diff(x);
y_s = diff(y);
err_s=sqrt(err(1:length(err)-1).^2 + err(2:length(err)).^2);

mv_graph(x,y, err, x_label, y_label,g_label)

mv_msg('Define area for background: click left and right limits');

%----- Find data inside of defined region
%  bkgd_x = all the data inside the background region

[x_cut,y_cut]=ginput(2);
bkgd_x=find(x > min(x_cut) & x < max(x_cut));
   
if isempty(bkgd_x);
   mv_msg('    No action taken: the whole graph is background!')
   return;
end
   
n=length(bkgd_x);

%  Integrated background is the mean of the counts inside the selected region
%  multiplied by the total width of the scan
bkgd = mean(y(bkgd_x))*(max(x)-min(x));
bkgd_err=mean(sqrt(err(bkgd_x).^2))*(max(x)-min(x));

%  Area of each trapezium is the height(of one point)*width
%  + 1/2*width*height difference
area = sum(x_s.*(y(1:length(y)-1) + 0.5*y_s));
area_error=sqrt(sum(x_s.*(err(1:length(y)-1).^2 + err_s.^2)));

area2 = num2str(area-bkgd);
area2_err=num2str(sqrt(area_error^2+bkgd_err^2));

msg = strvcat(['Area = ', num2str(area),' ± ', num2str(area_error), ', no background'],...
              ['Area-bkgd = ', area2, ' ± ', area2_err, ', Integ Bkgd = ', num2str(bkgd),' from ',num2str(n), ' pts']);
mv_msg(msg);

%----- Reset radio buttons

mv_rtidy(0)

return



