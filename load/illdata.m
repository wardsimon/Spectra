function [data,datastr,pres,pscan,stepvec,stepstr,com,header,head,flip,extrasearch]=illdata(filename,searchfield)
%
% function [data,datastr,pres,pscan,stepvec,stepstr,header,data_head,flip,extrasearch]=illdata(filename,searchfield)
%
% MATLAB function to read an ILL triple axis data file
%
% DFM 15.5.96 rev EF 07.07.97 (update)  rev ARW 08.09.98 (polarization)
%
% Reads ILL files (parameters and data) and extra fields if given

%--------- Initialize arrays ------------------------------

if nargin < 2
  searchfield = '';
end

data=[];datastr=[];pres=[];pscan=[];
stepvec=[];stepstr=[];com=[];header=[]; flip = 0;
head = []; extrasearch = [];

%--------- Open data file ---------------------------------

[fid,message]=fopen(filename,'r');
if (fid<0)
   fprintf(1,'ERROR on %s open: %s\n', filename, message);
   return;
end

%---- Read column headers ---------------------------------

fpos=ffind(filename,'DATA_:');
if isempty(fpos)
  %fprintf(1,'%s does not contain any data point.\n', filename);
  return
end
fseek(fid,fpos,'bof');
r=fgets(fid);
head=fgetl(fid);
head=['  ' head '  '];

%-----Read data -------------------------------------------
datapos = ftell(fid);
r=fgets(fid);
while length(r) > 2
   a=sscanf(r,'%f');
   if ~isempty(a)
   data=[data ; a']; end
   r=fgets(fid);
end

%----Process Header to extract labels ---------------------

datastr = [];
p = 1;
lh = length(head);
while p < lh
  letpos = find(isletter(head(p:lh)));
  if ~isempty(letpos)
    istart = letpos(1)+p-1;
    spapos = find(isspace(head(istart:lh)));
    if ~isempty(spapos)
      iend = spapos(1)+istart-2;
      toadd = head(istart:iend);
      if ~isempty(toadd) & find(isletter(toadd))
        if isempty(datastr)
          datastr = toadd;
        else
          datastr = str2mat(datastr,toadd);
        end
        p = iend+1;
      end
    else
      p = lh;
    end
  else
    p = lh;
  end
end

fpos=ffind(filename,'COMND:');
fseek(fid,fpos,'bof');
com=fgets(fid);

fseek(fid,0,-1); % rewind;
header = fread (fid,datapos-1);
header = setstr(header');


%----- Now read the resolution parameters

if nargout > 3

   pres=1:42;

   fpos=findstr(header,'PARAM');
   if isempty(fpos)
  disp('No param section.');
   else
  section=header(fpos(1):length(header));
    tosearch=str2mat('DM=','DA=',...
    'ETAM=','ETAA=','ETAS=',...
    'SM=','SS=','SA=',...
    'KFIX=','FX=',...
    'ALF1=','ALF2=','ALF3=','ALF4=',...
    'BET1=','BET2=','BET3=','BET4=',...
    'AS=','BS=','CS=',...
    ' AA=','BB=','CC=',...
    'AX=','AY=','AZ=',...
    'BX=','BY=','BZ=');
  [n,c]=size(tosearch);
  x=1:n;
    for j=1:n
    item=deblank(tosearch(j,:));
    fpos=findstr(section,item);
    if isempty(fpos)
      x(j) = NaN;
    else
      fpos=fpos(1)+length(item);
                        t = sscanf(section(fpos:length(section)),'%f');
      if isempty(t),
        disp([ 'WARNING: Could not read ' tosearch(j,:) ' value (invalid string)' ])
        t=0;
                        end
      x(j) = t;
    end
  end
  if x(3) == 0, x(3) = 30; disp('ETAM changed from 0 to 30 min'); end
  if x(4) == 0, x(4) = 30; disp('ETAS changed from 0 to 30 min'); end
  if x(5) == 0, x(5) = 30; disp('ETAA changed from 0 to 30 min'); end
  pres(1:30)=x;
   end

   fpos=findstr(header,'POSQE');
   if isempty(fpos)
  disp('No posqe section.');
   else
  section=header(fpos(1):length(header));
    tosearch=str2mat('QH=','QK=','QL=','EN=');
  [n,c]=size(tosearch);
  x=1:n;
    for j=1:n
    item=deblank(tosearch(j,:));
    fpos=findstr(section,item);
    if isempty(fpos)
      x(j) = NaN;
    else
      fpos=fpos(1)+length(item);
      x(j) = sscanf(section(fpos:length(section)),'%f');
    end
  end
  pres(31:34)=x;
   end

   x=[];
   fpos=ffind(filename,'STEPS');
   fseek(fid,fpos,'bof');
   r=fgets(fid);
% check if we have a second STEPS line
   r2 = fgets(fid);
   if ~isempty(strmatch('STEPS', r2))
    r = [ r, r2 ];
   end
   steps = [ ' ' r ' ' ];
   lsteps = length(steps);
   stepstr2 = '';

   dqh = 0;
   for id = findstr(steps,'QH');
  	 if isempty(id), break; end
     id1=findstr(steps((id(1)+3):lsteps),'=');
		 if isempty(id1), break; end
     if lsteps <= id(1)+id1(1)+3, break; end
     x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
     dqh = x(1);
     if isempty(stepstr2)
       stepstr2 = 'QH';
     else
       stepstr2 = str2mat(stepstr2,'QH');
     end
     stepvec = [ stepvec ; x(1) ];
   end

   dqk = 0;
   for id = findstr(steps,'QK');
     if isempty(id), break; end
     id1=findstr(steps((id(1)+3):lsteps),'=');
		 if isempty(id1), break; end
     if lsteps <= id(1)+id1(1)+3, break; end
     x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
     dqk = x(1);
     if isempty(stepstr2)
       stepstr2 = 'QK';
     else
       stepstr2 = str2mat(stepstr2,'QK');
     end
     stepvec = [ stepvec ; x(1) ];
   end

   dql = 0;
   for id = findstr(steps,'QL');
     if isempty(id), break; end
     id1=findstr(steps((id(1)+3):lsteps),'=');
		 if isempty(id1), break; end
     if lsteps <= id(1)+id1(1)+3, break; end
     x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
     dql = x(1);
     if isempty(stepstr2)
       stepstr2 = 'QL';
     else
       stepstr2 = str2mat(stepstr2,'QL');
     end
     stepvec = [ stepvec ; x(1) ];
   end

   den = 0;
   for id = findstr(steps,'EN');
     if isempty(id), break; end
     id1=findstr(steps((id(1)+3):lsteps),'=');
		 if isempty(id1), break; end
     if lsteps <= id(1)+id1(1)+3, break; end
     x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
     den = x(1);
     if isempty(stepstr2)
       stepstr2 = 'EN';
     else
       stepstr2 = str2mat(stepstr2,'EN');
     end
     stepvec = [ stepvec ; x(1) ];
   end


   pres(35:38)=[ dqh dqk dql den];
   pres(39:42)=[ 0 0 1 1];

   if isempty(stepvec)
     streq=findstr('=',r);
     for i=1:length(streq)
  stepvec =[stepvec sscanf(r(streq(i)+1:length(r)),'%f',1)];
     end
   end

%----Process STEPS header to extract labels ---------------------
   if ~isempty(stepstr2)
  stepstr = stepstr2;
   else
  p=findstr(steps,' ');
  K=p(filter([1 -1],1,p)>1); % start of non blank string analysing differences
  stepstr=[];
  for i=1:3:length(K)-1
    addstr = findstr(steps((K(i)+2):length(steps)),' ');
    addstr = min(addstr(1)+(K(i)+2),lsteps);
    addstr=steps((K(i)+1):addstr);
    if (i==1) stepstr=addstr;
    else
      stepstr=str2mat(stepstr,addstr);
    end
  end
   end
end

%----- Find if polarisation was used

ipol=findstr(header,'PAL');
if isempty(ipol)
  if ~isempty([ findstr('f1',lower(head)) findstr('f2',lower(head)) ])
    ipol = 1;
  end
end
if ~isempty(ipol)
  flip=1;
else
  flip=0;
end

%----- Calculate the scan limits

if nargout > 4

   pscan(1:4)=pres(31:34);
   [npts,dummy]=size(data);
   pscan(5:8)=pres(31:34)+pres(35:38)*(npts-1);

end

%----- Get extra optional fields

if ~isempty(searchfield)
  if iscell(searchfield)
    searchfield = str2mat(searchfield);
  end
  tosearch = searchfield;
  [n,c]=size(tosearch);
  x=1:n;
    for j=1:n
    item=deblank(tosearch(j,:));
    fpos=findstr(header,item);
    if isempty(fpos)
      x(j) = NaN;
    else
      fpos=fpos(1)+length(item);
      x(j) = sscanf(header(fpos:length(header)),'%f');
    end
  end
  extrasearch=x;
else
  extrasearch = [];
end

fclose(fid);
