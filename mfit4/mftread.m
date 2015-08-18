function [x, err]=mftread(file,par)
%
% function [x, err]=mftread(file,par)
% Read parameter with name par from .mft file name file

x=[];
err=[];
par=[' ' par ' '];

%------------- Open data file---------------------
fid=fopen(file,'r');				
if (fid<0) 
	error(['Couldn''t find ' file]);
end

line=fgetl(fid);

while line~=-1

	if strcmp(strtok(line,' '),'par')
		p=findstr(par,line);
		if ~isempty(p)
			t=sscanf(line,'%*s%*d%*s%f%f',2);
			x=[x; t(1)];
			err=[err; t(2)];
		end
	end
	line=fgetl(fid);
end

fclose(fid)

