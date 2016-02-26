function s_out = combine(toll,varargin)
%
% function s_out = combine(toll,s1,s2,....sn,'method','methods','indexing','relative')
%
% @SPEC1D/COMBINE function to combine two or more spectra.
%
% If the x values of two points differ by less
% than tolerance 'toll', then the points are combined.
%
% Depending on the optional 'method', points are combined as
% 'mean'    : Simple means for x and y, errors are averaged in quadrature.
% 'counts'	: Restablishes normalisation and original counts assuming
%             square-root statistics, is correct for normalised counts
% 'weight'	: Weights to inverse error. For more general data.
% Default is 'counts'
%
% Depending on the optional 'indexing', points are indexed as
% 'relative'    : Histogram bining of size 'toll'. This is the same as
%                 using the rebin function except for multiple spectum.
% 'absolute'	: The gap between the points is  always greater then 'toll'
%                 before any averaging.
% 'legacy'      : Replicate the original @spec1d/combine
% Default is 'relative' due to speed considerations.
%
% s1,s2,... can be single spectra or arrays of spectra.
%
% Example:
% Combine s1,s2 and s3 when x values differ by less than 0.01.
% >r = combine(0.01,s1,s2,s3) % combine s1,s2,s3 by counts method
% >r = combine(0.01,s1,s2,s3,'method','mean') % combine s1,s2,s3 by mean method
% >r = combine(0.01,s1,s2,s3,'method','weight','indexing',absolute) % combine s1,s2,s3 by
% weight method and indexing 'absolute'

% Simon Ward 25/01/2016


s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin(s_ind) = [];

p = inputParser;
p.CaseSensitive = false;
p.addParamValue('toll',@(x) isnumeric(x) && isreal(x));
p.addParamValue('s_in',@iscell);
p.addParamValue('indexing','relative',@ischar)
p.addParamValue('method','counts',@ischar);

p.parse(toll,s,varargin{:});
s = [p.Results.s_in{:}];

bin = p.Results.toll;
method = p.Results.method;
indexing = p.Results.indexing;

x = vertcat(s(:).x);
[x, ind ] = sort(x(:));
y = vertcat(s(:).y); y = y(ind);
e = vertcat(s(:).e); e = e(ind);
y_fit = vertcat(s(:).yfit); 
if ~isempty(y_fit)
    y_fit = y_fit(ind);
end

switch lower(indexing(1))
    case 'r'
        % Relative
        ind = [1; ceil(cumsum(diff(x(:)))/bin)]; % This is a faster way of [~,~,ind] = histcounts(x(:),min(x):bin:max(x));
    case 'a'
        % Absolute
        ind = zeros(size(x));
        ind(1) = 1;
        xtemp = x(1);
        for i = 2:length(x)
            if isempty(xtemp)
                xtemp = x(i);
            else
                xtemp = [xtemp x(i)]; %#ok<AGROW>
            end
            if diff(xtemp([1 end])) > bin
                ind(i) = ind(i-1)+1;
                xtemp = [];
            else
                ind(i) = ind(i-1);
            end
        end
    case 'l'
        % legacy
        s_out = combine_legacy(method,bin,s);
        return
    otherwise
        error('spec1d:combine:NotValidIndexing','%s is not a valid indexing method. See documentation',p.Results.method)
end

switch lower(method(1))
    case 'm'
        % Simple mean - Remastered to be GPU compatible
        N = accumarray(ind(:),ones(size(x)),[],@sum,NaN);
        xs = accumarray(ind(:),x(:),[],@sum,NaN)./N(:);
        ys = accumarray(ind(:),y(:),[],@sum)./N(:);
        es = sqrt(accumarray(ind(:),e(:).^2,[],@sum))./N(:);
        if ~isempty(y_fit)
            y_fit_s = accumarray(ind(:),y_fit(:),[],@sum)./N(:);
        else
            y_fit_s = [];
        end
    case 'c'
        % Counts
        % Reconstruct the montor count.
        % For zero error, this is impossible, so we take average monitor value.
        m = zeros(size(e),class(e)); % For GPU work.
        if any(y < 0)
            warning('spec1d:combine:NegativeY','Some Y values are negative. Using method "counts" is not recomended.')
            m(e~=0) = abs(y(e~=0))./e(e~=0).^2;
        else
            m(e~=0) = y(e~=0)./e(e~=0).^2;
        end
        if any(e==0)
            m(e==0) = mean(m(e~=0));
        end
        ms = accumarray(ind(:),m(:),[],@sum,NaN);
        ms(ms < eps) = eps; % Check for possibility of rounding errors
        xs = accumarray(ind(:),x(:).*m(:),[],@sum,NaN)./ms;
        ys = accumarray(ind(:),y(:).*m(:),[],@sum)./ms;
        if ~isempty(y_fit)
            y_fit_s = accumarray(ind(:),y_fit(:).*m(:),[],@sum)./ms;
        else
            y_fit_s = [];
        end
        es = (accumarray(ind(:),(e(:).*m(:)).^2,[],@sum).^0.5)./(ms.*accumarray(ind(:),ones(size(x)),[],@sum,NaN));
    case 'w'
        % Weight
        ms = accumarray(ind(:),1./e(:).^2,[],@sum);
        ms(ms < eps) = eps; % Check for possibility of rounding errors
        xs = accumarray(ind(:),x(:)./e(:).^2,[],@sum,NaN)./ms;
        ys = accumarray(ind(:),y(:)./e(:).^2,[],@sum)./ms;
        if ~isempty(y_fit)
            y_fit_s = accumarray(ind(:),y_fit(:)./e(:).^2,[],@sum)./ms;
        else
            y_fit_s = [];
        end
        es = 1./(sqrt(accumarray(ind(:),1./e(:).^2,[],@sum)).*accumarray(ind(:),ones(size(x)),[],@sum,NaN));
    otherwise
        error('spec1d:combine:NotValidMethod','%s is not a valid method. See documentation',p.Results.method)
end

% This takes into account the possibility of non-continuous indexing
r = s(1);
r.x = xs(~isnan(xs));
r.y = ys(~isnan(xs));
r.e = es(~isnan(xs));
r.yfit = y_fit_s;

s_out = spec1d(r);

