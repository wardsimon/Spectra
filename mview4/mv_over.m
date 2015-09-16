function mv_over
%
% MATLAB function to overlay buffers
% in mview
%
% DFM 21.4.97
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs < 2

   mv_msg(' Must select two or more buffers to overlay.');
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

%----- Perform operation on buffers: overlay

%----- All of data will be written as a large array of size
%      N x 4M, where N is the largest number of data points in
%      a scan and M is the number of scans. For a particular 
%      value of M, the columns are x, y, err, blank.

%----- Obtain y offset

def={num2str(0)};
offset=inputdlg({'Enter the offset in Y:'},'Overlay Offset',1,def);
if isempty(offset), return; end
offset=str2num(offset{1});

x=buffers(oper_buffs(1)).xobs;
y=buffers(oper_buffs(1)).yobs;
x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

N=length(buffers(1).xobs);

for i=2:noper_buffs

   x=buffers(oper_buffs(i)).xobs;
   y=buffers(oper_buffs(i)).yobs+offset*(i-1);
   x_min_new = min(x);
   x_max_new = max(x);
   y_min_new = min(y);
   y_max_new = max(y);

   if x_min_new <= x_min; x_min = x_min_new; end
   if x_max_new >= x_max; x_max = x_max_new; end
   if y_min_new <= y_min; y_min = y_min_new; end
   if y_max_new >= y_max; y_max = y_max_new; end
 
   y_max=y_max*1.1;
   N=max(N,length(buffers(1).xobs));  
 
end

%----- Determine N and M and create a matrix 

M=noper_buffs;

stored_scans=ones([N,4*M])*nan;

%----- Plot first data buffer

x  = buffers(oper_buffs(1)).xobs;
y  = buffers(oper_buffs(1)).yobs;
err= buffers(oper_buffs(1)).err;
x_label= buffers(oper_buffs(1)).x_label;
y_label= buffers(oper_buffs(1)).y_label;

mv_graph(x,y, err, x_label, y_label,'Overlay');

%----- Set new x and y ranges
hxy=get(findobj('Tag','mv_fixax'),'Checked');

if strcmp(hxy,'off')
   xrange=x_max-x_min;
   yrange=y_max-y_min;
   x_lim=[x_min-0.01*xrange x_max+0.01*xrange];
   y_lim=[y_min-0.1*yrange  y_max+0.1*yrange];
else
   x_lim=get(gca,'Xlim');
   y_lim=get(gca,'Ylim');
end

set(gca,'Xlim',x_lim); 
set(gca,'Ylim',y_lim);

[hmv_ctrl, hmv_data]=mv_figs;
figure(hmv_data);
hax=get(hmv_data,'CurrentAxes');  
axes(hax); delete(legend);
limits=[get(hax,'Xlim') get(hax,'Ylim')];

Marker_Size=str2num(get(findobj('Tag','mv_MarkerSize'),'String'));

%----- Loop over number of buffers

leg_h = zeros(1,noper_buffs);
leg_s = cell(1,noper_buffs);

for i=1:noper_buffs

   x  = buffers(oper_buffs(i)).xobs;  
   y  = buffers(oper_buffs(i)).yobs+offset*(i-1);
   err= buffers(oper_buffs(i)).err;

   stored_scans(1:length(x),1+(i-1)*4:i*4)=[x y err ones(size(x))];

%-------- Plot data from 2nd buffer onwards
 
   if i~=1
      mcolour={'c','y','g','b','k','r',...
               'c','y','g','b','k','r',...
               'c','y','g','b','k','r',...
               'c','y','g','b','k','r'};
      lso={'square','diamond','hexagram','x','pentagram','o',...
           'square','diamond','hexagram','x','pentagram','o',...
           'square','diamond','hexagram','x','pentagram','o',...
           'square','diamond','hexagram','x','pentagram','o'}; 
      hsel=line(x,y,'Marker',lso{i-1},...
                    'MarkerSize',Marker_Size,...
                    'Linestyle','None',...
                    'MarkerFaceColor',mcolour{i-1},...
                    'color',mcolour{i-1});

   	leg_h(i) = hsel;
	if ~isempty(buffers(oper_buffs(i)).datafile)
		leg_s{i} = buffers(oper_buffs(i)).datafile;
	else
		leg_s{i} = ['Buffer ' get(tmv_radio(i),'String') ];
	end

   end


%----- Write stored scans to userdata

end

set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',0);
set(hmv_data,'Userdata',stored_scans)
mv_uplot;

leg_h(1) = findobj('Tag','mv_selected');
leg_s{1} = buffers(oper_buffs(1)).datafile;

h = get(findobj('tag','mv_LegendStyle'),'string');
if isempty(h), h = 'on'; end
switch lower(h)
case {'on','best'}
    h = 0;
  case {'off','none'}
    h = [];
  case {'aside','outside','bestoutside'}
    h =-1;
  otherwise
    h = []
  end

if ~isempty(h), legend(leg_h, leg_s, h); end


%----- Reset radio buttons

mv_rtidy(0)

return

