function mv_func2(operation)
%
% MATLAB function to perform operation on two buffers.
% in mview.
%
% Allowed operations are: Divide, Interpolate, Multiply, Subtract
%
% DFM 21.4.97
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

x_length_test = length(buffers(oper_buffs(1)).xobs);

%----- Perform test on buffers to make sure they are the same length
%x_length_test = length(buffers(oper_buffs(1)).xobs);
%for i=2:noper_buffs
%   buff_index=oper_buffs(i);
%   x=buffers(buff_index).xobs;
%   x_length_test_new=length(x);
%   if x_length_test_new ~= x_length_test & ~strcmp(operation,'interpolate')
%      mv_msg(' Different x values: Interpolate first!')
%      return
%   end
%end

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

%----- Check that buffers overlap
if min(x1) > max(x2) | max(x1) < min(x2)
  mv_msg(' Scans don''t overlap');
  return
end

mflag=0;
if mean(m1)-mean(m2)~=0
  mv_msg(' WARNING: monitors in two scans are DIFFERENT!')
    mflag=1;
end
%----- Establish a tolerance for an operation
if length(x1) > 1
   toll=abs(x1(length(x1)) - ...
            x1(length(x1)-1))/4;
else
   toll=abs(x1(length(x1)))/10;
end

def={num2str(toll)};
tollin=inputdlg({'Input the tolerance (X axis):'},'Tolerance',1,def);
toll=str2num(tollin{1});

%----- Find values of x1, x2 that fit within tolerance
x1_pos=[]; x2_pos=[];
for i=1:length(x1)
  if find(abs(x1(i)-x2)<toll)
    x1_pos=[x1_pos;i];
    x2_pos=[x2_pos; find(abs(x1(i)-x2)<toll)];
  end
end

y=zeros([length(x_length_test),1]);
err=zeros([length(x_length_test),1]);
g_label=[];

%----- Perform one of the allowed operations
switch operation

   case 'divide'

   x1=x1(x1_pos); x2=x2(x2_pos);
   y1=y1(x1_pos); y2=y2(x2_pos);
   e1=e1(x1_pos); e2=e2(x2_pos);
   m1=m1(x1_pos); m2=m2(x2_pos);
      x   = x1;
      izero=find(y2==0);
      if ~isempty(izero)
         mv_msg('Removing zero valued data in denominantor');
         x(izero) =[];
         y1(izero)=[]; y2(izero)=[];
         e1(izero)=[]; e2(izero)=[];
      end
      y   = y1./y2;
      err = sqrt((y1.*e2).^2+(y2.*e1).^2)./y2.^2;
    if mflag, m1=ones(size(y)); end;
    mon=m1;
      x_label= x_label_1;
      y_label= 'Counts';
      g_label= ['(' g_label_1 ')' ' / ' '(' g_label_2 ')'];
      label='DIV: ';

   case 'interpolate'

      if length(x1)==length(x2)
         if isempty(find(x1~=x2))
            mv_msg(' No interpolation:  X ranges are identical');
            return
         end
      end
      in=1;
      while min(x2(in)) < min(x1); x2(1)=[]; in=in+1; end
      in=1;
      while max(x2(in)) > max(x1); x2(end)=[]; in=in+1; end
      x     = x2;
      y     = interp1(x1,y1,x2);
      err   = interp1(x1,e1,x2);
    mon   = interp1(x1,m1,x2);
      x_label = x_label_1;
      y_label = 'Counts';
      g_label = ['Interpolated: ' '(' g_label_1 ')' ' on ' '(' g_label_2 ')'];
      label='INTER: ';

   case 'multiply'

    x1=x1(x1_pos); x2=x2(x2_pos);
    y1=y1(x1_pos); y2=y2(x2_pos);
    e1=e1(x1_pos); e2=e2(x2_pos);
    m1=m1(x1_pos); m2=m2(x2_pos);
      x   = x1;
      y   = y1 .* y2;
    if mflag, m1=ones(size(y)); end;
    mon=m1;
      err = sqrt((y1.*e2).^2+(y2.*e1).^2);
      x_label= x_label_1;
      y_label= 'Counts';
      g_label= ['(' g_label_1 ')' ' * ' '(' g_label_2 ')'];
      label='MULTI: ';

   case 'subtract'

    x1=x1(x1_pos); x2=x2(x2_pos);
    y1=y1(x1_pos); y2=y2(x2_pos);
    e1=e1(x1_pos); e2=e2(x2_pos);
    m1=m1(x1_pos); m2=m2(x2_pos);
      x    = x1;
      y    = y1-y2;
      err  = sqrt((e1).^2+(e2).^2);
    if mflag, m1=ones(size(y)); end;
    mon=m1;
      x_label= x_label_1;
      y_label= 'Counts';
      g_label= ['(' g_label_1 ')'  ' - ' '(' g_label_2 ')'];
      label='SUB: ';

end

%----- Update buffers

[incr_buffs]=mv_ubuff(x,y,err,mon,x_label,y_label,g_label,[]);
set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',incr_buffs);

%----- Change text in window to file name and plot results

hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
for i=1:noper_buffs; label=[label, num2str(oper_buffs(i)) '  ']; end
set(hmv_text(incr_buffs),'String',label,'ToolTipString',label);

%----- Plot results

mv_graph(x,y,err,x_label,y_label,g_label);
mv_msg(label);

%----- Reset radio buttons

mv_rtidy(0)

return
