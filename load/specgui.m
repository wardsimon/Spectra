function [x,y,err,xlab,ylab,monitor]=specgui(filename)
%
%function [x,y,err,xlab,ylab,monitor]=specgui(filename)
%
%   This function extracts scans from data files of the type generated at
%   the SPEC. Search scan variable and ask user for 
%   the X,Y variables to import. EF 08.07.97 DFM 12.6.96
%
%Output variable
%   x        - the independent variable
%   y        - the dependent variable
%   err      - uncertainty in dependent variable = sqrt(y)
%   xlab     - name of 'x' data
%   ylab     - name of 'y' data
%
% This routine either accepts a simple filename input 
% argument, or alternatively a compound file
% specification of the form:
%
%      filename,{options,...}
%
% Valid options are =
%      X=string     name for X, use '-' for auto setting
%      Y=string     name for Y
%      M=string     name for monitor, can also be : none, use '-' for auto setting
%      E=string     name for Y error, can also be : sqrt(y) or none
%      N=number     normalisation value
%      S=number     scan number or field numbers ('*' for all, '-' for main field)
%      gui          ask user with GUI
%      setpar       set parameters if possible (in rescal or mfit)
%     
% example :
%      batch('foo.txt,X=2:10,M=-,E=sqrt(y),Y=Counts,N=norm(y)')

% uses : ffind.c as a mex file. specatch.m

%----- Parse filespec --------------------------------------
filespec = filename;
filename2 = filename;
[fspec filespec]=strtok(filespec,',');
while ~isempty(filespec)
   [s filespec]=strtok(filespec,',');
   fspec=str2mat(fspec,s);
end
[nargs,nchars]=size(fspec);

i=strmatch('S=',fspec);
if ~isempty(i)
	scan=deblank(fspec(i(1),3:nchars)); 
else
	scan = '';
	hs=findobj('Tag','mfload_scan');
	if ~isempty(hs) scan=get(hs,'String'); end
	if isempty(scan) 
		scan = 1; 
		disp('Warning : SPEC load : What scan should I load in file ?');
		disp('Warning :             Going to scan #1 (option ''S=1'')...');
		scan = 1;
	end
end

[x, y, err, xlab, ylab,monitor]=specbatch([ filename2 ',gui' ]);
disp('TIP : select "File/Reload" to reach a new scan in SPEC file.');
