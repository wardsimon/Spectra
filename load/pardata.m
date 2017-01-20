function [data,datastr]=pardata(file)
%
% function pardata(file)
% 
% MATLAB function to read parameters from mfit parameter file
%
% DFM 27.4.97

%------------- Open data file---------------------


data=[];
datastr=[];

fid=fopen(file,'r');				
if (fid<0) 
	return
end

text=[' Reading parameters from file ' file];
fprintf('\n %s \n\n',text)

first=1;
chisq2=[];	
while (1)
       
   line=fgetl(fid);
   if ~isempty(line) & (line==-1), break, end
   
   x=[];
   while strcmp(strtok(line,' '),'par')		
      t=sscanf(line,'%*s%*d%*s%f%f');      
      x=[x,t'];
      
      if first==1		
         datastr=str2mat(datastr,sscanf(line,'%*s%*d%s%*f%*f'));	
	 datastr=str2mat(datastr,['err_' sscanf(line,'%*s%*d%s%*f%*f')]);		
      end 

      line=fgetl(fid);
   end
   if strcmp(strtok(line,'Chi^2'),'% Norm. ')
      chisq2=[chisq2; sscanf(line,'%*s%*s%*s%f')];      
   end
   data=[data; x];
   first=isempty(data);	

end
data=[data, [chisq2 sqrt(chisq2)]];
datastr(1,:)=[];
datastr=str2mat(datastr,'Chisq^2','err_Chisq^2');
fclose(fid);



