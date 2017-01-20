function rc_getp(option,infile)
%
% MATLAB function to get a parameter file for RESCAL.
%
% Note: If option == 'pfile', update parameters in main window
%       If option == 'ifile', update parameters in instrument parameter window
%       If option == 'sfile', update parameters in simulation parameter window
%       If option == 'tfile', update parameters in scan part of trixpar window
%
% DFM 9.5.95
%

if strcmp(option,'pfile')

   ftype='*.par';

elseif strcmp(option,'sfile')

   ftype='*.sim';

elseif strcmp(option,'ifile')

   ftype='*.cfg';

end

%----- Either prompt for file name or use infile
file=[];

if nargin == 1 

   [file,fdir]=uigetfile(ftype,'Rescal parameter file');

   str=[fdir file];

else 

   str = infile;
   
end

%if file~=0

   fid=fopen(str,'r');				
   if (fid<0) 
      datafile=0;
      return
   end

%end

%-------------- Initialize arrays-----------------

data=[];		
header='';
text=fgetl(fid);

%----------------- Load data-----------------------

while (text>0)
   [temp count]=sscanf(text,'%f');			
   if (isempty(temp)) 
      header=[header text];
   else
      if (count==size(data,2) | isempty(data))
         data=[data; temp'];
      end
   end
      text=fgetl(fid);
end
fclose(fid);

p=data(:,1);

%----- Write parameters to edit windows

if strcmp(option,'pfile')

%----- Parameters for main RESCAL window

         if length(p) ~=42
            disp(' Warning: Incorrect number of rescal parameters')
         end 

         hrc_paras=get(findobj('Tag','hrc_paras'),'Userdata');
         hrc_units_paras=get(findobj('Tag','hrc_units_paras'),'Userdata');
         hrc_units_current=get(findobj('Tag','hrc_units_current'),'Userdata');
	 
         for i=1:length(hrc_paras)
            set(hrc_paras(i),'String',num2str(p(i)));
         end

	 % also set the energy units to Ang-1
	 set(hrc_units_paras(1),'Value',1);
	 set(hrc_units_current(1),'Value',1);

elseif strcmp(option,'ifile')

%----- Parameters for instrument window

         if length(p) ~=27
            disp(' Warning: Incorrect number of rescal parameters')
         end 

         hrc_inst_paras=get(findobj('Tag','hrc_inst_paras'),'Userdata');

         set(hrc_inst_paras(1),'Value',p(1)+1);
         set(hrc_inst_paras(2),'String',num2str(p(2)));
         set(hrc_inst_paras(3),'String',num2str(p(3)));
         set(hrc_inst_paras(4),'Value',p(4)+1);
         set(hrc_inst_paras(5),'String',num2str(p(5)));
         set(hrc_inst_paras(6),'String',num2str(p(6)));
         set(hrc_inst_paras(7),'Value',p(7)+1);
         set(hrc_inst_paras(8),'String',num2str(p(8)));
         set(hrc_inst_paras(9),'String',num2str(p(9)));
         set(hrc_inst_paras(10),'String',num2str(p(10)));
         set(hrc_inst_paras(11),'Value',p(11)+1);
         set(hrc_inst_paras(12),'String',num2str(p(12)));
         set(hrc_inst_paras(13),'String',num2str(p(13)));

         for i=14:length(hrc_inst_paras)
            set(hrc_inst_paras(i),'String',num2str(p(i)));
         end

%----- Update options in instrument window

         rc_inst('Source');
         rc_inst('Guide');
         rc_inst('Sample');
         rc_inst('Detector');

elseif strcmp(option,'sfile')

%----- Parameters for Monte-Carlo window

         hmc_paras=get(findobj('Tag','hmc_paras'),'Userdata');

         for i=1:length(hmc_paras)
            set(hmc_paras(i),'String',num2str(p(i)));
         end

elseif strcmp(option,'tfile')

         htrix_scan_paras=get(findobj('Tag','htrix_scan'),'Userdata');

         for i=1:length(htrix_scan_paras)
            set(htrix_scan_paras(i),'String',num2str(p(i)));
         end

   end





