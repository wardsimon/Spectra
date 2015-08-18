function [data,datastr]=ntasdata(file)
%
% function [data,datastr]=ntasdata(file)
%
% MATLAB function to load NEW TASCOM file, returning
% the whole  file as an array 'data', with headings 'datastr'
%
% DFM: 27.4.97
%
%----- Open TASCOM data file

data=[]; datastr=[];

dferror=0;

fid=fopen(file,'r');

if (fid <0) 

   datafile=0;
   dferror=1;
   return;
   
end  

%----- Read through data file 

test='zzzz';
data=[];
datastr=[];

while strncmp(test,'#fip',4)==0

    dataline=fgetl(fid);
    if ~isempty(dataline)
	if dataline==-1, dferror=1; return, end 
    end
    test=strtok(dataline);
    if ~isempty(test), test=dataline(1:4); end
           
end

%----- Create column headers

[kill,dataline]=strtok(dataline);
while ~isempty(dataline)

   [column_lab,dataline]=strtok(dataline);
   datastr=strvcat(datastr,column_lab);

end

%----- Read data

while(1)

   dataline=fgetl(fid);
   if dataline==-1, 
      return 
   elseif ~isempty(dataline)
      a=sscanf(dataline,'%f'); data=[data; a'];
   end

end

fclose(fid)

return
