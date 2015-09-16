function [x, y, err, xlab, ylab,monitor]=mcabatch(filespec)
%
%function [x, y, err, xlab, ylab,monitor]=mcabatch(filespec)
%
% MVIEW/MFIT load routine for MCA text data files. 
% This routine either accepts a simple filename input 
% argument, or alternatively a compound file
% specification of the form:
%
% filename,{options,...}
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


[x, y, err, xlab, ylab,monitor]=multibatch([ filespec ',Y=2:ncolumns,X=1' ]);
if isempty(x)
	return
end
n = length(y); % last number is sweep number

xlab = 'Channels';
ylab = 'Counts';

if n
	i = log(n)/log(2);
	if (i == round(i))
		mon = y(n);
		y = y(1:(n-1));
		x = x(1:(n-1));
		ylab = [ ylab sprintf(' (%i sweeps)',mon) ];
	else
		mon = 1;
	end
end

% transfert to mfit params
hmf_pars = findobj('Tag','mf_ParWindow');
h=[];
if ~isempty(hmf_pars)
	h=get(hmf_pars,'userdata');
end
[npars d]=size(h);
p=1:npars;
for i=1:npars
	c=str2num(get(h(i,1),'string'));
	if (isempty(c))
		p(i) = 0;
	else
		p(i) = c;
	end
	pnames = lower(get(h(i,3),'string'));
	if ~isempty(findstr(pnames,'sweeps')) | ~isempty(findstr(pnames,'scans')) | ~isempty(findstr(pnames,'monit'))
		p(i) = mon;
	end
end
% make modifs

for i=1:npars
	set(h(i,1),'string',num2str(p(i),6));
end
