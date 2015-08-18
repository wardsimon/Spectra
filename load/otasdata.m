function [data,datastr]=otasdata(file)
%
% function [data,datastr]=otasdata(file)
%
% MATLAB function to load OLD TASCOM file, returning
% the whole  file as an array 'data', with headings 'datastr'
%
% DFM: 27.4.97
%
%----- Open TASCOM data file

data = []; datastr=[];

dferror=0;

fid=fopen(file,'r'); 

if (fid<0)

   datafile=0;
   dferror=1;
   return;

end

%----- Determine the number of data columns
   
nparas  =fscanf(fid, '%g %g', [2 1]);
nskip   =nparas(1);
ncolumns=nparas(2);

%----- Read headers

datastr=[];
i=0;
while i < nskip-1

   dataline=strtok(fgetl(fid)); 
   if ~isempty(dataline)
      if i~=0
         datastr=strvcat(datastr,dataline);
      end
      i=i+1;
   end
      
end 

datastr(1:nskip-ncolumns-2,:)=[];
   
%----- Read the names of the column headers   

length_test=fscanf(fid,'%g',[1 inf]);
ndata_points=length(length_test)/ncolumns;
data=reshape(length_test,ncolumns,ndata_points)';

fclose(fid);

return
