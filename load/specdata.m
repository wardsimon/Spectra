function [data, datastr, com, head]=specdata(filename, scan)
%
%function [data, datastr, com, head]=specdata(file, scan)
%
% This function enables to load SPEC scan data
%
% MZ 8.3.95
% DFM 30.9.95
% EF 4.9.97

data=[]; datastr=''; com = ''; head = '';

if nargin < 2, scan=[]; end
if isempty(scan)
  scan=1;
  disp('Warning : SPECdata : no scan number precised : using scan #1');
end

%---- Read through data file to find right scan number (uses ffind mex file)

fpos=ffind(filename,['#S ' num2str(scan) ' ']);
if fpos<0
   errordlg(['Scan ' num2str(scan) ' not found'],'Spec data load error:')
   return
end

fid=fopen(filename);
if isempty(fpos) fpos=1; end
if fpos<0 fpos=0; end
fseek(fid,fpos,'bof');
com=fgets(fid);                   % This line contains the scan command issued

%---- Read column headers ---------------------------------
t='zz';
while strcmp(t,'#L') == 0
   r=fgetl(fid);
   if (length(r) > 1) t=r(1:2); end
end
head=['  ' r(4:max(size(r))) '  '];

%-----Read data -------------------------------------------
data=[];
r=fgets(fid);
while (max(size(r))>2 & r(1)~='#')
   a=sscanf(r,'%f');
   data=[data ; a'];
   r=fgets(fid);
end
fclose(fid);

%----- Write column headers to a matrix -----------------------

K=findstr(head,'  ');
datastr=[];
for i=1:length(K)-1
   datastr=strvcat(datastr, fliplr(deblank(fliplr(deblank(head(1,K(i)+1:K(i+1)-1))))));
end





