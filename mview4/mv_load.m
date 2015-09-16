function [x,y,err,xlab,ylab]=mf_load(file)
%
% MFIT function [x,y,err]=mf_load(file)
%	MZ 29.11.94
%
% This is a basic load routine for MFIT, illustrating the required
% syntax. The routine takes the name of a data file (including path) as a
% parameter and returns the column vectors x, y, err.
% On error the routine must exit with datafile=0.

%------------- Open data file---------------------
fid=fopen(file,'r');				
if (fid<0) 
	datafile=0;
	return
end

%-------------- Initialize arrays-----------------
data=[];		
header='';
text=fgetl(fid);

%----------------- Load data-----------------------
while (text>0)
	[temp count]=sscanf(text,'%f');			
	if (temp==[]) 
		header=[header text];
	else
		if (count==size(data,2) | data==[])
			data=[data; temp'];
		end
	end
	text=fgetl(fid);
end
fclose(fid);

x=data(:,1);
y=data(:,2);
[dummy,ncols]=size(data)
if ncols > 3
   err=data(:,3);
else
   err=sqrt(data(:,2));
end
xlab='x';
ylab='y';

end
