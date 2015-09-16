function [x,y,err,xlab,ylab,mon]=multi(file)
%
% MFIT load routine [x,y,err,xlab,ylab,monitor]=multi(file)
%
% This is the basic MFIT multicolumn load routine. When called it presents
% a window allowing the user to choose which columns to plot. The data file
% must be a white-spaced multicolumn ascii file, with optional header text.
% If the last line of the header is a list of column headings, these will
% appear in the menu list.
%
% MZ 5.6.95, DFM 2.5.97

%===== Parse file spec ============================================
x=[];y=[];err=[];xlab='';ylab='';

[file s]=strtok(file,',');
[xcol s]=strtok(s,',');
[ycol s]=strtok(s,','); 
[ecol s]=strtok(s,',');
[mcol s]=strtok(s,',');

if ~isempty(xcol) & ~isempty(ycol)
   xcol=str2num(xcol);
   ycol=str2num(ycol);
   if isempty(ecol)
      ecol='sqrt(y)';
   elseif ~isempty(str2num(ecol))
      ecol=str2num(ecol);
   end
   if isempty(mcol)
      mcol='none';
   elseif ~isempty(str2num(mcol))
      mcol=str2num(mcol);
   end
end

%============= Load data ===========================================

%----- Open file ---------------------------------------------------

disp(['Loading ' file]);
fid=fopen(file,'r');				
if (fid<0) 
   disp('file not found');
   return; 
end;
header='dummy';
text=fgetl(fid);
istext=1;

%------Read header and first line of data --------------------------

while (istext==1)                            
   [row1 count]=sscanf(text,'%f');	% See if is numbers... 
   if (isempty(row1))                             % ...if not...
      header=str2mat(header, text);          % ...add string to header.
      text=fgetl(fid);
   else
      istext=0;
   end
end
ncol=length(row1);                           % number of columns in data file

%------ Read rest of data and reshape into matrix ------------------

data=fscanf(fid,'%f');                       % Vectorized load of rest of data
nrow=length(data)/ncol;                      % Work out number of rows
if (nrow*ncol~=length(data))                 
   error('Bad data format');
end
data=[row1'; reshape(data,ncol,nrow)'];      % Reshape vector to matrix
fclose(fid);                                 % close input file

%------- Choose columns to load -------------------------------------

collist=header(size(header,1),:);            % Extract last line of header
datastr=[];
while ~isempty(collist)

   [column_lab,collist]=strtok(collist);
   datastr=strvcat(datastr,column_lab);

end
collist=datastr;

if isempty(xcol)
   [xcol, ycol, ecol, mcol, xlab, ylab]=...  % Pass it to column choosing dialog routine
   mf_coldg(collist,ncol,[],[],[],[],1,1,[], [ filename ' Column Selector' ]);       % with actual number of columns
else
   xlab='x';
   ylab='y';
end

%========= Extract relevant data columns ==================================	
if ~isempty(xcol)

%---- x and y data --------------------------------------------------------

   x=data(:,xcol);                              % x data 
   y=data(:,ycol);                              % y data

%---- work out y uncertainty ----------------------------------------------

   if (ecol==ncol+1)
   	err=sqrt(abs(y));                              % ...root(y) errors requested        
   elseif (ecol>=ncol+2)
      err=1e-3*max(abs(y))*ones(size(y));       % no error bars - so equal weights
   else
   	err=data(:,ecol);                         % specified errors
   end

%---- normalize to monitor?-------------------------------------------------

   if mcol <= ncol                      % Assume counting statistics 
      mon=data(:,mcol);                         % in monitor
      err=(y./mon).*sqrt(1./(abs(y)+1)./mon);
      y=y./mon;
      disp('Warning: assumed counting statistics in monitor');
   else
      mon = ones(size(y));
   end
end
