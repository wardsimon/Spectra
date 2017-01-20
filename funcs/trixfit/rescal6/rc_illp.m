function rc_illp

%
% MATLAB function to get a parameter file from ILL data format.
%
% DFM 9.5.95
%

   ftype='*.scn';

%----- prompt for file name

   [file,fdir]=uigetfile(ftype,'Rescal parameter file');
str=[fdir,file]

if file~=0

   fid=fopen(str,'r');				

end

%----------------- Load data-----------------------


% Read parameters from ILL files.
p=[1:42];
x=rc_lineread(str,'PARAM');
i=[1:1:30];
j=[1,2,16,17,24,3,4,5,7,6,8,9,10,11,12,13,14,15,18,19,20,21,22,23,25,26,...
27,28,29,30];
p(i)=x(j);
x=rc_lineread(str,'POSQE');
p(31:34)=x;
x=rc_lineread(str,'STEPS');
p(35:38)=x;
p(39:42)=[ 0 0 1 1];

%----- Parameters for main RESCAL window

         if length(p) ~=42
            disp(' Warning: Incorrect number of rescal parameters')
         end 

         hrc_paras=get(findobj('Tag','hrc_paras'),'Userdata');

         for i=1:length(hrc_paras)
            set(hrc_paras(i),'String',num2str(p(i)));
         end

end

