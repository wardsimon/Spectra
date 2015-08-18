function tomfit(x,y,err,selected,p,fixed)
% MFIT data send
% function tomfit(x,y,err,selected,p,fixed);
% 		Send data to MFIT
% 		EF 15.11.97
%
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
if nargin == 1
	y=x;
	x=[];
end
if nargin < 3, err = []; end
if nargin < 4, selected = []; end
if nargin < 5, p = []; end
if nargin < 6, fixed = []; end

update = 0;

if isempty(hmf_ctrl) | ~hmf_ctrl
	mfit;
end

if (isempty(hmf_data) | ~hmf_data) & ~isempty(y)
	mf_dwin;
	update = 1;
end

if (isempty(hmf_pars) | ~hmf_pars) & ~isempty(p)
	mf_newfn;
	update = 1;
end

pnames = [];
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

%---------- Get userdata in MFit for empty params ------------------


if ~isempty(hmf_data) & hmf_data
	userdata =get(hmf_data,'userdata');
	if ~isempty(userdata)
		if isempty(x) x=userdata(:,1); end
		if isempty(y) y=userdata(:,2); end
		if isempty(err) err=userdata(:,3); end
		if isempty(selected) selected=userdata(:,4); end
	end
end

if ~isempty(hmf_pars) & hmf_pars
	[i,j,k,l] = mf_rpars;
	if isempty(p) p=i; end
	if isempty(fixed) fixed = k; end
	pnames = l;
end

% now send data

if ~isempty(p) | ~isempty(fixed) 
	mf_upars(p,[],fixed);
end
if ~isempty(p) | ~isempty(pnames)
	mf_pwin([],pnames,p);
end

i = length(y);

if length(x) ~= i x=1:i; end
if length(err) ~= i, err = y/1000; end;
if length(selected) ~= i, selected = ones(i,1); end

x=x(:); y=y(:); err=err(:); selected = selected(:);
userdata = [ x y err selected ];

set(hmf_data,'userdata',userdata);
mf_gdata('noload');
set(hmf_data,'userdata',userdata);
if update
	mf_uplot('all')
else
	mf_uplot('sel')
end

