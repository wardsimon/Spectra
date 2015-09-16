function [p,pmat,x,y] = mf_autogs(x, y, selected, toset, pmat, pnamemat,funfile)
% MFIT auto guess routine
% [p,pmat,x,y] = mf_autogs(x, y, selected, toset, pmat,pnamemat,fcnname)
% determines start parameters for peaks fitting
% in MFIT signal, and retuns start params
% This routine searches for : background, peak position, width and intensity.
% toset is 1 (set) or 0 (don't set) vector for [ position intensity width ]
% pmat is [ index left_width right_width  max_pos Intensity Width ]

% uses : mf_fndpks, diff, mf_smooth, mf_rmspks


if nargin < 2
	x=[]; y=[];
end
if  nargin < 3
	selected = [];
end
if nargin < 4
	toset = [];
end

if nargin < 5
	pmat = [];
end

if nargin < 6
	pnamemat = [];
end

p=[]; pmat = [];

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
if (isempty(x) | isempty(y)) & ~isempty(hmf_data) & hmf_data
	data=get(hmf_data,'userdata');
	selected = data(:,4);
	x=data(:,1);
	y=data(:,2);
	if sum(selected) <= 1
		disp('Can''t guess : no points selected in Data Window');
		return
	else
		idx = find(~isnan(y) & ~isinf(y) & (selected ~= 0));
	end
end
if isempty(y)
	disp('Can''t guess : empty data !!');
	return;
end
if isempty(x)
	x=1:length(y);
end
if (length(selected)<length(x))
	selected = ones(1,length(x));
end
if (isempty(toset))
	toset = [ 1 1 1 ];
end

% search peaks in data -------------------------------

if isempty(pmat)
	pmat = mf_fndpks(x,y); % make y positive non zero
	% pmat = [ index left_width right_width  max_pos Intensity Width ]

	[dummy, pmat] = mf_rmspks(x,y,pmat);
end

% sort peaks by intensity -------------------------------

[p, i] = sort(-pmat(:,5));	% sort in descending order
pmat = pmat(i,:);
i = find(selected(pmat(:,1))); % get peaks in selected region
pmat = pmat(i,:);
iint = 1; iwidth = 1; ipos = 1; ibackg = 1;

% modify parameters -------------------------------

h=[];
if nargin<7 %~exist('funfile')
	funfile = findobj('Tag','mf_FitFuncFile');
	if ~isempty(funfile) 
		funfile = get(funfile,'string');
		disp([ 'Auto guess for function : ' funfile ]);
	else 
		funfile ='';
	end
end
if ~isempty(hmf_pars) & hmf_pars & isempty(pnamemat)
	h=get(hmf_pars,'userdata');
	[oldp,dp,fixed] = mf_rpars;
elseif isempty(pnamemat)
	disp('No parameter window opened...');
	p = [];
	return
end

if ~isempty(h)
	[npars d]=size(h);
else
	[npars d]=size(pnamemat);
	fixed = zeros(npars,1);
end
[nmat d] = size(pmat);
parmat = ones(npars,2);
for i=1:npars
	iall = [ iint ipos ];
	if ~isempty(h)
		c=str2num(get(h(i,1),'string'));
	else
		c = [];
	end
	if (isempty(c))
		p = 0;
	else
		p = c;
	end
	if ~isempty(h)
		pnames = lower(get(h(i,3),'string'));
		pfix = get(h(i,3),'value');
	else
		pnames = lower(pnamemat(i,:));
		pfix = 0;
	end
	n = [ findstr(pnames,'backg') findstr(pnames,'bkg') findstr(pnames,'const')];
	if (~isempty(n) & ~pfix)
		if (ibackg)
			p = min(y);
			if p == 0
				p = min(y(find(y>p)));
			end
			ibackg = 0;
			pfix = 0;
		else
			p = 0;
			pfix = 1;
		end
	end
	n = [ findstr(pnames,'energ') findstr(pnames,'pos') findstr(pnames,'cent') findstr(pnames,'freq')  ];
	if ~isempty(n) & (ipos <= nmat)
		if (max(ipos - iall) > 1) % then need to stop peak setting
			break;
		else
			if (toset(1) & ~pfix) p = pmat(ipos,4); end
			if p > min(x) & p < max(x)
				ipos = ipos +1;
			end
		end
	end
	n = [ findstr(pnames,'wid') findstr(pnames,'damp') findstr(pnames,'gamma')];
	if ~isempty(n) & (iwidth <= nmat)
		if (max(iwidth - iall) <= 1) % then peak damping setting
			if (toset(3) & ~pfix) p = pmat(iwidth,6)/2; end
			iwidth = iwidth +1;
		end
	end
	n = [ findstr(pnames,'int') findstr(pnames,'heigh') findstr(pnames,'amp') ];
	if ~isempty(n) & (iint <= nmat)
		if (max(iint - iall) > 1) % then need to stop peak setting
			break;
		else
			n = find(pmat(iint,4) <= x);
			if isempty(n)
				n = pmat(iint,1);
			end
			n = n(1);
			if (toset(2) & ~pfix) p = y(n); end
			iint = iint +1;
		end
	end
	parmat(i,:) = [ p pfix ];
end
% make modifs
if isempty(pnamemat)
	i = find(fixed);
	if (~isempty(i))
		parmat(i,1) = oldp(i)';
	end
	fprintf(1,'Auto Guessed about %i peaks\n',max([ iint iwidth ipos]) - 1);
	mf_upars(parmat(:,1),[],parmat(:,2))
	mf_uplot('fit');
	p = parmat;
else
	p = parmat(:,1); % pfix = 0
end

