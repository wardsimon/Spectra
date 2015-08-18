function mv_save(save_comm)
%
% MATLAB function to save buffer to an ascii file
%
% DFM 27.4.97
%

%----- Unload data from invisible storage

buffers=get(findobj('Tag','hmv_buffers'),'Userdata');
[dummy,nbuffs]=size(buffers);

%----- Unload handles of radio buttons

tmv_radio=get(findobj('Tag','tmv_radio'),'Userdata');
noper_buffs=length(tmv_radio);

if noper_buffs~=1

   mv_msg(' Must select only one buffer to save.');
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

%----- Perform operation on buffers: save

buff_sv = oper_buffs;

cur_dir=cd;

output_data=[buffers(buff_sv).xobs buffers(buff_sv).yobs buffers(buff_sv).err];

%----- Get current output directory

cur_out_dir=get(findobj('Tag','mv_OutDir'),'String');
cur_out_file=get(findobj('Tag','mv_OutFile'),'String');

format long e
if nargin==0
   [out_file, out_dir]=uiputfile([cur_out_dir '*.dat'],'Save buffer as');
   if out_file ==0; return; end
   set(findobj('Tag','mv_OutDir'),'String',out_dir);
   set(findobj('Tag','mv_OutFile'),'String',out_file);
else
   out_dir=cur_out_dir;
   out_file=cur_out_file;
end
%eval(['save ' [out_dir out_file] ' output_data -ascii'])
%save([out_dir out_file],'output_data','-ascii')
t=fix(clock);
fid = fopen([out_dir out_file],'w');
fprintf(fid,'#Mview data file saved %d.%d.%d   %d:%d:%d\n',t(3:-1:1),t(4:6));
fprintf(fid,'#datafile: %s\n' ,buffers(buff_sv).datafile);
fprintf(fid,'#title  : %s\n' ,buffers(buff_sv).g_label);
fprintf(fid,'#x_label: %s\n', buffers(buff_sv).x_label);
fprintf(fid,'#y_label: %s\n', buffers(buff_sv).y_label);
for i=1:size(output_data,1)
    fprintf(fid,'%G ',output_data(i,:));
    fprintf(fid,'\n');
end

fclose(fid)

mv_msg(['Saving to file: ' out_dir out_file ])
format

%----- Reset radio buttons

mv_rtidy(0)

return

