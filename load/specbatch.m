function [x, y, err, xlab, ylab,monitor]=specbatch(filespec)
%
%function [x, y, err, xlab, ylab,monitor]=specbatch(file)
%
% MVIEW/MFIT load routine for SPEC data files. 
% This routine either accepts a simple filename input 
% argument, or alternatively a compound file
% specification of the form:
%
%      filename,{options,...}
%

% Valid options are =
%      X=string     name for X, use '-' for auto setting, '#' for point number
%      Y=string     name for Y
%      M=string     name for monitor, can also be : none, use '-' for auto setting
%      E=string     name for Y error, can also be : sqrt(y) or none
%      N=number     normalisation value
%      S=number     scan number
%      gui          ask user with GUI
%      setpar       set parameters if possible (in rescal or mfit)
%     
% example :
%      batch('foo.txt,X=2,M=-,E=sqrt(y),Y=Counts,N=norm(y),S=5,gui')     


% MZ 11.10.95 and DFM 2.11.95, EF 4.09.97
%
% this file is composed of those different parts :
% * filespec analysis, retrieval of options
% * load data and datastr (plus possible options)
% * display file info for user                 <- file format dependent
% * extract column guess or choice              <- file format dependent
% * possibly call GUI                          <- programer's choice
% * last column check
% * extract x,y,error,monitor,labels
% * options : automatic set params,...         <- file format dependent
%
% for GUI, you MIGHT use mf_coldg (x,y,error,monitor selector) at your choice.


%===== Parse filespec, open data file, and get column names ====================

%---- Set default column names and flags

scan='';
xname='none';
yname='Detector';
mname='none';
ename='-';
normf=1;

setpars = 0;	% automatic set of rescal pars
gui = 0;		% no gui choice

x=[]; y= []; err=[]; ylab=''; xlab=''; monitor = [];

%----- Parse filespec --------------------------------------

[fspec filespec]=strtok(filespec,',');
while ~isempty(filespec)
   [s filespec]=strtok(filespec,',');
   fspec=str2mat(fspec,s);
end
[nargs,nchars]=size(fspec);

%----- Update scan parameters from filespec---------------------------

i=strmatch('X=',fspec);
if ~isempty(i)
	xname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('Y=',fspec);
if ~isempty(i)
	yname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('M=',fspec);
if ~isempty(i)
	mname=deblank(fspec(i(end),3:nchars)); end

i=strmatch('N=',fspec);
if ~isempty(i)
	normf=deblank(fspec(i(end),3:nchars)); end

i=strmatch('E=',fspec);
if ~isempty(i)
	ename=deblank(fspec(i(end),3:nchars)); end

i=strmatch('S=',fspec);
if ~isempty(i)
	scan=deblank(fspec(i(end),3:nchars)); end

filename=deblank(fspec(1,:));

if strmatch('gui',fspec)
	gui = 1;
end

if strmatch('setpar',fspec)
	setpars=1;
end

scan = str2num(scan);

%===== Get data from file ============================

disp([ 'Loading SPEC data : ' filename ]);
[data,datastr,com,head]=specdata(filename, scan);
if isempty(data) return; end		% load error
if isempty(datastr)
	[nchar,ncolumns] = size(data);
	datastr=[];
	for i=1:ncolumns
		datastr=str2mat(datastr,sprintf('col%i',i));
	end
	datastr(1,:)=[];
end
[ncolumns,nchar] = size(datastr);
ndatastr = '';
for i=1:ncolumns  				% set a header line
	ndatastr = [ ndatastr ' ' datastr(i,:) ];
end

% ===== print informations ============================

disp(com);
p=findstr(head,'  ');                  % posns of column names in head
%---- Check for temperature stability ---------------------
pp=findstr(head,'  DegK  ');
T=[];
if ~isempty(pp)
	tcol=find(pp==p);
	T=mean(data(:,tcol));
	sT=std(data(:,tcol));
	fprintf(1,'T=%2f +/- %2f K\n',T,sT)
end

[n,c] = size(data);
fprintf(1,'Data is ( %ix%i ) in scan #%s.\n',n,c,num2str(scan));
fprintf(1,'%s\n',ndatastr);
nc = min(10,ncolumns);
for i=1:min(2,n)
	fprintf(1,'%g ',data(i,1:nc));
	fprintf(1,'\n');
end
fprintf(1,'    ...\n');
for i=max(1,n-1):n
	fprintf(1,'%g ',data(i,1:nc));
	fprintf(1,'\n');
end

%===== extract x,y,mon columns guesses ============================

n = ncolumns;
xcol=eval(xname,'strmatch(xname,datastr,''exact'')');
%if isempty(xcol) xcol = ncolumns+1; end

if (strcmp(yname,'-'))
	yname='Detector';
end
ycol=eval(yname,'strmatch(yname,datastr,''exact'')');

if (strcmp(mname,'-'))
	mname = 'Monitor';
end
mcol=eval(mname,'strmatch(mname,datastr,''exact'')');
%if isempty(mcol) mcol = ncolumns+1; end

if (strcmp(ename,'-'))
	ename='sqrt(y)';
end
ecol=eval(ename,'strmatch(ename,datastr,''exact'')');
if isempty(ecol)
	if strcmp(ename,'sqrt(y)') ecol = ncolumns+1; end
end


%===== GUI option call  ============================
if gui
	[xcol, ycol, ecol, mcol, xlab, ylab, scantmp, normf,keep]=mf_coldg(datastr, ncolumns,xcol,ycol,ecol,mcol,scan,normf, [], [ filename ' Column Selector' ]);
	if isempty(xcol) return; end 	% cancel button
	if isempty(scantmp) scantmp='1'; end
	if (~strcmp(scantmp,scan)) % need to reach a new part of file 
		scan = scantmp;
		if (xcol<=ncolumns) 
			xcol = datastr(xcol,:); 
		else
			xcol='#';
		end
		if (ycol<=ncolumns) ycol = datastr(ycol,:); 
		else ycol = num2str(ycol); end
		if (mcol<=ncolumns) mcol = datastr(mcol,:); 
		else mcol = 'none'; end
		if (ecol<=ncolumns) ecol = datastr(ecol,:); 
		else
			if (ecol == ncolumns+1) ecol = 'sqrt(y)';
			else ecol = 'none'; end
		end
		filespec = sprintf('%s,X=%s,Y=%s,M=%s,E=%s,S=%s,N=%s',filename,xcol,ycol,mcol,ecol,scan,normf);
		if setpars
			filespec = [ filespec ',setpars' ];
		end
		if (keep)
			filespec = [ filespec ',gui' ];
		end
		[x, y, err, xlab, ylab]=specbatch(filespec);
		return
	end

end

%===== last check for columns =============

if isempty(ecol) ecol=ncolumns+2; end
if isempty(mcol) mcol=ncolumns+1; end

if (isempty(xcol)) %----- Find first column that varies for x
	xcol=find(std(data));
	xcol=xcol(1);
end

if (isempty(ycol)) %----- Find first column that varies for y next to x
	ycol=find(std(data(:,(xcol+1):ncolumns)));
	ycol=ycol(1);
end
%===== extract x,y values ============================

%---- Work out which columns to extract, and extract --------

%----- Select columns to analyse

% first group columns in case x/ycol is a group of columns

if (xcol <= ncolumns)
	x=data(:,xcol);
	[n,c]=size(x);
	x=reshape(x',n*c,1);
	ncx=n*c;
else
	ncx = -1; % causes x = 1:length(y)
end

y=data(:,ycol);
[n,c]=size(y);
y=reshape(y',n*c,1);

if (ncx ~= n*c)
	x=1:(n*c);
	x=x';
end


if mcol <= ncolumns
	monitor=data(:,mcol);
else
	monitor=ones(length(x),1);
end
[n,c]=size(monitor);
monitor=reshape(monitor',n*c,1);

[monzeros]=find(monitor==0);	%----- Test to see if selected monitor is zero
if ~isempty(monzeros);
	disp(' ')
	disp('  Warning: selected monitor has some zeros')
	disp(' ')
	monitor=ones(length(x),1);
end

if (ecol == ncolumns+1)
	err=sqrt(abs(y));			% ...root(y) errors requested        
elseif (ecol >= ncolumns+2)
	err=1e-3*max(abs(y));		% no error bars - so equal weights
else
	err=data(:,ecol);			% specified errors
end
[n,c]=size(err);
err=reshape(err',n*c,1);
err = err./monitor;

%----- Normalise counts

normfv=eval(num2str(normf));

y=normfv*y ./monitor; 
err=normfv*err;

%----- Create labels

   if (xcol <= ncolumns)
   	xlab=datastr(xcol(1),:); 
   else
	xlab='Point number';
   end
   xlab=sprintf('%s (scan #%d)',xlab,scan);
   if ~isempty(T)
	xlab=sprintf('%s T=%.2f K',xlab,T);
   end
   ymon=[' ( per  ' num2str(mean(monitor)) '  Monitor ) * ' num2str(normf)];
   ylab2=ymon;
   ylab1=fliplr(deblank(fliplr(deblank(datastr(ycol(1),:)))));
   ylab=[ylab1 ylab2];

%===== Set params ============================

if setpars

% no set par yet

end % if setpars
