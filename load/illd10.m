function [x,y,err,xlab,ylab,monitor]=illd10(filename)
%
% function [x,y,err,xlab,ylab,monitor]=illd10(filename)
%
% MATLAB function to read an ILL d10 data file
%
% DFM 14.5.96
%

%--------- Initialize arrays ------------------------------
x=[];y=[];err=[];xlab='';ylab='';
scan_info=[]; data=[];

%--------- Open data file ---------------------------------

fid=fopen(filename,'r');
if (fid<0)
   error('File not found');
end

%---- Find scan info

fpos=ffind(filename,'valco:');
fseek(fid,fpos,'bof');
dummy=fgets(fid);

%----- Read Scan info

r=fgets(fid);
while max(size(r)) > 2
   a=sscanf(r,'%f');
   scan_info=[scan_info ; a'];
   r=fgets(fid);
end

%---- Find data 

fpos=ffind(filename,'Numor');
fseek(fid,fpos,'bof');
dummy=fgets(fid);

%----- Read Scan info

r=fgets(fid);
while max(size(r)) > 2
   a=sscanf(r,'%f');
   data=[data ; a'];
   r=fgets(fid);
end

%----- Select data

y=data(:,3);
err=sqrt(abs(y));

qstart=[scan_info(1,1) scan_info(1,2) scan_info(1,3)];
qend=[scan_info(5,2) scan_info(5,3) scan_info(5,4)];
deltaq=qend-qstart;

qvary=find(deltaq);
qstep=deltaq(qvary)/(length(y)-1);
qvec=qstart(qvary):qstep:qend(qvary);

x=qvec';

%----- Create labels

ylab='Counts';
labels=str2mat('H','K','L');
xlab=labels(qvary);
monitor = ones(size(y));
