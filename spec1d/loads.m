function s=loads(varargin)
%
% function s=loads(filetype,datafile)
%
% MATLAB function to load datafile of type filetype  
% to a spec1d object s. Here varargin contains the 
% string used to load the required single file (scan),  
% or multiple files (scans).
%
% Notes:
% 1. Filetype must be a valid routine from the mfit/mview load library 
%    These include:  
%
%    xyeload     :  x,y,error
%    multibatch  :  multicolumn
%    specbacth   :  SPEC file
%    tasbatch    :  TASCOM file
%    illbatch    :  ILL 3-axis file
%    parbatch    :  MFIT parameter output file
% 
% 2. After first use it is only necessary to specify the part of the
%    input string that changes. If variables are to be removed at a 
%    later time, then set them equal to '[]', e.g.  '...,M=[],...'
%    will remove the monitor.
%
% 3. Loads multiple files, and multiple scans (SPEC data files)
%       [n1 n2 n3] - load files or scans specified by n1, n2 and n3 into a spec1d array
%       [n1:n3]    - same as above
%       [n1:n2:n3] - load files or scans from n1 to n2 with a regular spacing of n2
%       [n1;n3]    - combine(n1,n3) into single spec1d object
%       [n1++n3]   - combine(n1,n2,n3) into single spec1d object
%
% Examples:  
% 1. Load the 3rd  and then 4th scan  from SPEC file cesb.05
%    >>s3=loads('specbatch','cesb.05,X=K,Y=Detector,M=Monitor,S=3');    
%    >>s4=loads('S=4')
% 
% 2. Load the TASCOM file cuge108.dat
%    >>s=loads('tasbatch','cuge108.dat,X=OM,Y=I,M=Mon');
%
% 3. Create an array of spectra, e.g. [108, 109, 110, 111, combine(112,114), combine(115,116,117)]
%    >>s=loads('tasbatch','cuge[108 109:111 112;114 115++117].dat,X=OM,Y=I,M=Mon'); 
%
% 4. Load scans 20 through 24 into spec1d array
%    >>s=loads('specbatch','upd3.05,X=Theta,Y=Exp_Hutch,M=Monitor,S=[20:24]'); 
%
% 5. Load scans 46,51,56,....161,166
%    >>s=loads('specbatch','cdreo7.02,X=Exp_Hutch,Y=DegK,M=[],S=[46:5:166]');
%
% Version 2.0, April 2001
% Des McMorrow and Henrik Ronnow

%----- Define globals for storing file details

global S1D_filetype
global S1D_datafile
global S1D_datafilestring
global S1D_readVars
global S1D_readVals

%----- filetype and datafile details specified

if nargin==2
    
   filetype=varargin{1};
   if ~exist(filetype); disp('Error: filetype not valid'); return; end
   S1D_filetype=filetype;
   datafilestring=varargin{2};
   S1D_datafilestring=datafilestring;
   
%----- datafile details only; check if filetype was given in earlier call   
   
elseif nargin==1
    
   if isempty(S1D_filetype)
      disp('Error: must specify file type to load')
      return
   else
       filetype=S1D_filetype;
       datafilestring=varargin{1};
   end
    
end    

%----- parse datafile details 

datafields='';
rem=datafilestring;
while (any(rem))   
   [stripstring,rem]=strtok(rem,',');  
   datafields=strvcat(datafields,stripstring);   
end
   
temp=datafields;
temp(any((temp=='=')'),:)=[];
if isempty(temp) & isempty(S1D_datafile)
   disp('Error: must specify file name properly')
   return
elseif  ~isempty(temp)  
   S1D_datafile=deblank(temp);
end   

%----- parse datafile columns to read, and check which ones have changed

readVars=[];
readVals=[];
[l1,l2]=size(datafields);
for i=1:l1
   if ~isempty(findstr('=',datafields(i,:)))
   [c1,c2]=strtok(datafields(i,:),'=');
   readVars=strvcat(readVars,deblank(c1));  
   readVals=strvcat(readVals,deblank(c2(2:end)));   
   end    
end

if isempty(S1D_readVars)
   S1D_readVars=readVars;
   S1D_readVals=readVals;
else
   New_readVars=S1D_readVars;
   New_readVals=S1D_readVals;
   
%----- Change vals of specified vars and append new ones  

[l1,l2]=size(readVals);
   for i=1:l1
      im=strmatch(readVars(i,:),New_readVars);
      if ~isempty(im)
         a=deblank(readVals(i,1:end));
         [l3,l4]=size(New_readVals);
         New_readVals(im,1:l4)=blanks(l4);
         New_readVals(im,1:length(a))=a;
      else 
         New_readVals=strvcat(New_readVals,deblank(readVals(i,1:end)));
         New_readVars=strvcat(New_readVars,deblank(readVars(i,1:end)));   
      end   
   end    
   
   S1D_readVars=New_readVars;
   S1D_readVals=New_readVals;
   
end
%---- Kill vars with a val of []   

[l1,l2]=size(S1D_readVals);
i=1;
while i <= l1   
   im=strmatch('[]',S1D_readVals(i,:));
   if ~isempty(im)
      S1D_readVars(i,:)=[];
      S1D_readVals(i,:)=[];
   else    
      i=i+1;
   end
   [l1,l2]=size(S1D_readVals);
end    

%----- Construct datafilestring

datafilestring=deblank(S1D_datafile);
readVars=S1D_readVars;
readVals=S1D_readVals;
[l1,l2]=size(readVars);
for i=1:l1
   datafilestring=[datafilestring ',' deblank(readVars(i,:)) '=' deblank(readVals(i,:))];
end

%----- Check to see if multiple data files or scans are to be loaded using [...] notation

multiscans=0;
iss=strmatch('S',S1D_readVars);
if ~isempty(iss)
 if strmatch(S1D_readVals(iss),'['); multiscans=1; end   
end   

if ~isempty(findstr(S1D_datafile,'['))

   datafile=S1D_datafile;
   specs=datafile(findstr(datafile,'[')+1:findstr(datafile,']')-1);

   [evalstring,expnumbers]=specsparse(specs);

   filefront=datafilestring(1:findstr(datafile,'[')-1);
   fileback=datafilestring(findstr(datafile,']')+1:end);

   ndigit=size(num2str(expnumbers(:)),2);
   for n=1:length(expnumbers)
      [x,y,err,xlab,ylab]=feval(filetype,[filefront,num2str(expnumbers(n),['%0.',num2str(ndigit),'i']),fileback]);
      if ~isempty(x)
      datafile=strtok(datafilestring,',');
      [x,nsort]=sort(x);
      ss.x=x;
      ss.y=y(nsort);
      ss.e=err(nsort);
      ss.x_label=xlab;
      ss.y_label=ylab;
      ss.datafile=datafile;
      ss.yfit=[];
      s(n)=spec1d(ss);
      end
   end
   eval(['s=',evalstring]);

%----- Check to see if multiple scans are to be loaded using [...] notation
      
elseif multiscans

   scans=S1D_readVals(strmatch('S',S1D_readVars),:);
   specs=scans(findstr(scans,'[')+1:findstr(scans,']')-1);

   [evalstring,expnumbers]=specsparse(specs);    

   filefront=datafilestring(1:findstr(datafilestring,'[')-1);
   fileback=datafilestring(findstr(datafilestring,']')+1:end);
   
   ndigit=size(num2str(expnumbers(:)),2);
   for n=1:length(expnumbers)
      [x,y,err,xlab,ylab]=feval(filetype,[filefront,num2str(expnumbers(n),['%0.',num2str(ndigit),'i']),fileback]);
      if ~isempty(x)
      datafile=strtok(datafilestring,',');
      [x,nsort]=sort(x);
      ss.x=x;
      ss.y=y(nsort);
      ss.e=err(nsort);
      ss.x_label=xlab;
      ss.y_label=ylab;
      ss.datafile=datafile;
      ss.yfit=[];
      s(n)=spec1d(ss);
      end
   end
   eval(['s=',evalstring]);   
   
%----- Load single file
        
else
   [x,y,err,xlab,ylab]=feval(filetype,datafilestring);
   if isempty(x); return; end
   datafile=strtok(datafilestring,',');
   [x,nsort]=sort(x);
   ss.x=x;
   ss.y=y(nsort);
   ss.e=err(nsort);
   ss.x_label=xlab;
   ss.y_label=ylab;
   ss.datafile=datafile;
   ss.yfit=[];
   s=spec1d(ss);
   
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evalstring,expnumbers]=specsparse(specs)
%
% function [evalstring,expnumbers]=specsparse(specs) 
%
% MATLAB function to parse specs information
%
% Henrik Ronnow, February 2001

rawnumbers={''};
rawoperators={''};
lastspec='';
for n=1:length(specs)
   if (double(specs(n))>47 & double(specs(n))<58) % its a number
       rawnumbers{end}=[rawnumbers{end} specs(n)];
      if strcmp(lastspec,'o') 
         rawoperators{end+1}='';
      end
      lastspec='n';
   else
      rawoperators{end}=[rawoperators{end} specs(n)];
      if lastspec=='n'
         rawnumbers{end+1}='';
      end
      lastspec='o';
   end
end
rawoperators=rawoperators(1:(length(rawnumbers)-1));

% Expand operators
for n=1:length(rawnumbers)
   numbers(n)=str2num(rawnumbers{n});
end
expnumbers=[];
operators={};
n=1;
while n<=length(rawoperators)
   if strcmp(rawoperators{n},':')
      if n<length(rawoperators)
        if strcmp(rawoperators{n+1},':')
          for m=numbers(n):numbers(n+1):(numbers(n+2)-numbers(n+1));
             expnumbers(end+1)=m;
             operators{end+1}=' ';
          end
          n=n+1;
        else
          for m=numbers(n):(numbers(n+1)-1);
             expnumbers(end+1)=m;
             operators{end+1}=' ';
          end
        end
     else
       for m=numbers(n):(numbers(n+1)-1);
          expnumbers(end+1)=m;
          operators{end+1}=' ';
       end
     end
   elseif strcmp(rawoperators{n},'++')
      if n<length(rawoperators)
        if strcmp(rawoperators{n+1},'++')
          for m=numbers(n):numbers(n+1):(numbers(n+2)-numbers(n+1));
             expnumbers(end+1)=m;
             operators{end+1}=';';
          end
          n=n+1;
        else
          for m=numbers(n):(numbers(n+1)-1);
             expnumbers(end+1)=m;
             operators{end+1}=';';
          end
        end
     else
       for m=numbers(n):(numbers(n+1)-1);
          expnumbers(end+1)=m;
          operators{end+1}=';';
       end
      end
   else
      expnumbers(end+1)=numbers(n);
      operators{end+1}=rawoperators{n};    
   end
   n=n+1;
end
expnumbers(end+1)=numbers(end);
operators{end+1}=' ';

evalstring='[';
n=1;
if operators{n}==';';
   evalstring=[evalstring 'combine([s(' num2str(n) ') '];
   incombine=1;
else
   evalstring=[evalstring 's(' num2str(n) ') '];
   incombine=0;
end
   
for n=2:length(expnumbers)
   if (incombine & operators{n}==';') | (~incombine & operators{n}~=';')
      evalstring=[evalstring ' s(',num2str(n),') '];
   elseif ~incombine & operators{n}==';'
      evalstring=[evalstring 'combine([s(' num2str(n) ') '];
      incombine=1;
   elseif incombine & operators{n}~=';'
      evalstring=[evalstring ' s(' num2str(n) ')]) '];
      incombine=0;
   end
end
evalstring=[evalstring '];'];
