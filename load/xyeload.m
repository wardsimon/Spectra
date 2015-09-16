function [x,y,err,xlab,ylab,monitor]=xyeload(file)
%
% MFIT function [x,y,err,xlab,ylab,monitor]=xyeload(file)
%	MZ 29.11.94
%
% This is a basic load routine for MFIT, illustrating the required
% syntax. The routine takes the name of a data file (including path) as a
% parameter and returns the column vectors x, y, err.
% On error the routine must exit with datafile=0.

%------------- Open data file---------------------
x=[]; y= []; err=[]; ylab=''; xlab='';

fid=fopen(file,'r');				
if (fid<0) 
	datafile=0;
	return
end

%-------------- Initialize arrays-----------------
data=[];		
header='dummy';
text=fgetl(fid);

%============= Load data ===========================================

%------Read header and first row of data --------------------------
istext=1;
while (istext==1)                            
	[row1 count]=sscanf(text,'%f');			
 	if isempty(row1) 
 		header=str2mat(header, text);
   	text=fgetl(fid);
 	else
		istext=0;
 	end
end

%------ Read data and reshape into matrix --------------------------
data=fscanf(fid,'%f');                       % Read data into vector (for speed)
ncol=length(row1);                           % Use row1 to work out how to reshape
nrow=length(data)/ncol;
if (nrow*ncol~=length(data))                 
	error('Bad data format');
end
data=[row1'; reshape(data,ncol,nrow)'];      % Reshape vector to matrix
fclose(fid);                                 % close input file

x=data(:,1);
y=data(:,2);
if length(data(1,:))==2
    err=sqrt(data(:,2));
else
    err=data(:,3);
end

%===== Make x and y column labels ===================================

%-----Try to be clever and make labels from last line of header ----------
s=deblank(header(size(header,1),:));                        % last line of header
i=find(s==9);                                      % replace tabs by spaces
s(i)=32*ones(size(i));
[xlab s]=strtok(setstr(s),' ');                    % x label is first word
[ylab s]=strtok(setstr(s),' ');                    % y label is second word
[elab s]=strtok(setstr(s),' ');                    % elabel is third word 

%----- Wrong number of columns? Then just 'x' and 'y'------------------
if ~(length(elab)>0 & length(s)==0)                % s should be empty if 3 words
   xlab='x';
   ylab='y';
end

monitor = ones(size(y));


