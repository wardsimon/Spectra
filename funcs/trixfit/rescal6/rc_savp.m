function [p]=rc_savp(option)
%
% MATLAB function rc_savp saves parameter  from RESCAL windows
%        either to a file or a
% Note: If option == 'pfile' , save RESCAL parameters to file *.par
%       If option == 'ifile' , save instrument parameters to file *.cfg
%       If option == 'sfile' , save simulation parameters to file *.sim
%       If option == 'spars' ,   unload simulation parametes to p
%       If option == 'respars' , unload all rescal parameters to p
%       If option == 'tpars'   , unload trix scan parameters to p
%
% DFM 9.5.95
%
% Updated so that it takes account of choice of energy units.
% Alan Tennant 10th April 1998

cur_dir=cd;
psim=[];ptrix=[];pres=[];pinst=[];
if strcmp(option,'pfile') | strcmp(option,'respars')

   ftype='*.par';

   pres=zeros([42,1]);

   hrc_paras=get(findobj('Tag','hrc_paras'),'Userdata');
   hrc_units_paras=get(findobj('Tag','hrc_units_paras'),'Userdata');

   for i=1:length(hrc_paras)
       pres(i)=str2num(get(hrc_paras(i),'String'));
   end

   % set energy to Angs-1
   kfix=pres(9);
   value=get(hrc_units_paras(1),'Value');

   % convert from current units to Angs-1
   if value==2  % meV->Angs-1
     K_angs=sqrt(kfix/2.072);
   elseif value==3  % THz->Angs-1
     K_angs=sqrt(kfix*4.1357/2.072);
   else                % Angs-1
     K_angs=kfix;
   end
   pres(9)=K_angs; % now the energy is in Angs-1

end

if strcmp(option,'ifile') | strcmp(option,'respars')

   ftype='*.cfg';

   pinst=zeros([27,1]);

   hrc_inst_paras=get(findobj('Tag','hrc_inst_paras'),'Userdata');

   pinst(1) =get(hrc_inst_paras(1),'Value');
   pinst(2) =str2num(get(hrc_inst_paras(2),'String'));
   pinst(3) =str2num(get(hrc_inst_paras(3),'String'));
   if pinst(1)  == 1, pinst(3)=pinst(2); end              % Circular source

   pinst(4) =get(hrc_inst_paras(4),'Value');
   pinst(5) =str2num(get(hrc_inst_paras(5),'String'));
   pinst(6) =str2num(get(hrc_inst_paras(6),'String'));
   if pinst(4)  == 1, pinst(5)=0; pinst(6)=0; end         % No guide

   pinst(7) =get(hrc_inst_paras(7),'Value');
   pinst(8) =str2num(get(hrc_inst_paras(8),'String'));
   pinst(9) =str2num(get(hrc_inst_paras(9),'String'));
   pinst(10)=str2num(get(hrc_inst_paras(10),'String'));
   if pinst(7)  == 1, pinst(10)=pinst(9); pinst(9)=pinst(8); end  % Cylindrical sample

   pinst(11)=get(hrc_inst_paras(11),'Value');
   pinst(12) =str2num(get(hrc_inst_paras(12),'String'));
   pinst(13) =str2num(get(hrc_inst_paras(13),'String'));
   if pinst(11) == 1, pinst(13)=pinst(12); end            % Circular detector

   pinst(1)=pinst(1)-1;   pinst(4)=pinst(4)-1;   pinst(7)=pinst(7)-1;   pinst(11)=pinst(11)-1;


   for i=14:length(hrc_inst_paras)
       pinst(i)=str2num(get(hrc_inst_paras(i),'String'));
   end
   
end


if strcmp(option,'sfile') | strcmp(option,'spars')

   ftype='*.sim';

   hmc_paras=get(findobj('Tag','hmc_paras'),'Userdata');

   for i=1:length(hmc_paras)

      psim(i)=str2num(get(hmc_paras(i),'String'));

   end

end

if strcmp(option,'tpars')

   htrix_scan_paras=get(findobj('Tag','htrix_scan'),'Userdata');

   for i=1:length(htrix_scan_paras)

      ptrix(i)=str2num(get(htrix_scan_paras(i),'String'));

   end

end

p=[pres; pinst; psim; ptrix];

if strcmp(option,'pfile') | strcmp(option,'sfile')

   [out_file, out_dir]=uiputfile(ftype,'Save parameters as');
   if out_file ==0; return; end

   eval(['cd ' out_dir]);

   fid=fopen(out_file,'w');				

   if (fid<0) 
      datafile=0;
      return
   end

   fprintf(fid,'%8.3f \n',p);
   fclose(fid);

   eval(['cd ' cur_dir]);

end 

if strcmp(option,'ifile')

   [out_file, out_dir]=uiputfile(ftype,'Save parameters as');
   if out_file ==0; return; end

   eval(['cd ' out_dir]);

   fid=fopen(out_file,'w');				

   if (fid<0) 
      datafile=0;
      return
   end

   sout=             '   % =0 for circular source, =1 for rectangular source'; 
   sout=str2mat(sout,'   % width/diameter of the source (cm)');
   sout=str2mat(sout,'   % height/diameter of the source (cm)');
   sout=str2mat(sout,'   % =0 No Guide, =1 for Guide');
   sout=str2mat(sout,'   % horizontal guide divergence (minutes/Angs)');
   sout=str2mat(sout,'   % vertical guide divergence (minutes/Angs)'); 
   sout=str2mat(sout,'   % =0 for cylindrical sample, =1 for cuboid sample');
   sout=str2mat(sout,'   % sample width/diameter perp. to Q (cm)');
   sout=str2mat(sout,'   % sample width/diameter along Q (cm)');
   sout=str2mat(sout,'   % sample height (cm)');
   sout=str2mat(sout,'   % =0 for circular detector, =1 for rectangular detector');
   sout=str2mat(sout,'   % width/diameter of the detector (cm)');
   sout=str2mat(sout,'   % height/diameter of the detector (cm)');
   sout=str2mat(sout,'   % thickness of monochromator (cm)');
   sout=str2mat(sout,'   % width of monochromator (cm)');
   sout=str2mat(sout,'   % height of monochromator (cm)');
   sout=str2mat(sout,'   % thickness of analyser (cm)');
   sout=str2mat(sout,'   % width of analyser (cm)');
   sout=str2mat(sout,'   % height of analyser (cm)');
   sout=str2mat(sout,'   % distance between source and monochromator (cm)');
   sout=str2mat(sout,'   % distance between monochromator and sample (cm)');
   sout=str2mat(sout,'   % distance between sample and analyser (cm)');
   sout=str2mat(sout,'   % distance between analyser and detector (cm)');
   sout=str2mat(sout,'   % horizontal curvature of monochromator 1/radius (m-1)');
   sout=str2mat(sout,'   % vertical curvature of monochromator (m-1)');
   sout=str2mat(sout,'   % horizontal curvature of analyser (m-1)');
   sout=str2mat(sout,'   % vertical curvature of analyser (m-1)');

   for i=1:length(p)
      fprintf(fid,'%8.3f %45s \n',p(i),sout(i,:));
   end

   fclose(fid);

   eval(['cd ' cur_dir]);

end 

return

