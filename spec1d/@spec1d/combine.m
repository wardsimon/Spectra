function s_out = combine(toll,varargin)
%
% function s_out = combine(toll,s1,s2,....sn,'method','methods','indexing','auto')
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
<<<<<<< HEAD
% 'relative'    : Histogram bining of size 'toll'. This is the same as
%                 using the rebin function except for multiple spectum.
%                 The end bin might not be a 'full' bin.
% 'histogram'   : Histogram bining of size 'toll'. This is the same as
%                 using the rebin function except for multiple spectum.
%                 We try to make the first and last element 'toll'/2.
% 'absolute'	: The gap between the points is  always greater then 'toll'
%                 before any averaging.
=======
%                   'auto'   The default 'auto' algorithm chooses a bin
%                            width to cover the data range and reveal the
%                            shape of the underlying distribution.
%                  'scott'   Scott's rule is optimal if X is close to being
%                            normally distributed, but is also appropriate
%                            for most other distributions. It uses a bin width
%                            of 3.5*STD(X(:))*NUMEL(X)^(-1/3).
%                     'fd'   The Freedman-Diaconis rule is less sensitive to
%                            outliers in the data, and may be more suitable
%                            for data with heavy-tailed distributions. It
%                            uses a bin width of 2*IQR(X(:))*NUMEL(X)^(-1/3),
%                            where IQR is the interquartile range.
%                'sturges'   Sturges' rule is a simple rule that is popular
%                            due to its simplicity. It chooses the number of
%                            bins to be CEIL(1 + LOG2(NUMEL(X))).
%                   'sqrt'   The Square Root rule is another simple rule
%                            widely used in other software packages. It
%                            chooses the number of bins to be
%                            CEIL(SQRT(NUMEL(X))).
>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8
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

% Simon Ward 19/01/2017


s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin(s_ind) = [];

p = inputParser;
p.CaseSensitive = false;
p.addRequired('toll',@(x) isnumeric(x) && isreal(x));
p.addRequired('s_in',@iscell);
p.addParameter('indexing','auto',@ischar)
p.addParameter('method','counts',@ischar);

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

<<<<<<< HEAD
switch lower(indexing(1))
    case 'r'
        % Relative
        if x(1)==x(2)
            x(2) = x(2)+eps;
        end
        ind = [1; ceil(cumsum(diff(x(:)))/bin)]; % This is a faster way of [~,~,ind] = histcounts(x(:),min(x):bin:max(x));
    case 'h'
        % Histogram
        ind = 1 + (ceil(x(end)/bin)*bin - floor(x(1)/bin)*bin)/bin - sum(bsxfun(@rdivide,x(:),linspace(floor(x(1)/bin)*bin-bin/2, ceil(x(end)/bin)*bin+bin/2, (ceil(x(end)/bin)*bin - floor(x(1)/bin)*bin)/bin+1))<1,2) +1;
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
=======
switch lower(indexing(1:2))
    case 'au'
        if isnan(bin)
            bin = fminsearch(@optim_bin,(max(x)-min(x))/max(ceil(log2(numel(x))+1),1));
>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8
        end
        edges = calc_bins(bin);
    case 'sc'
        edges = scottsrule(x, min(x), max(x), false);
    case 'fd'
        edges = fdrule(x, min(x), max(x), false);
    case 'st'
        edges = sturgesrule(x, min(x), max(x), false);
    case 'sq'
        edges = sqrtrule(x, min(x), max(x), false);
    case 'le'
        % legacy
        s_out = combine_legacy(method,bin,s);
        return
    otherwise
        error('spec1d:combine:NotValidIndexing','%s is not a valid indexing method. See documentation',p.Results.method)
end

<<<<<<< HEAD
=======
ind = arrayfun(@(i) i*((x >= edges(i)) & (x < edges(i+1))),1:(length(edges)-1),'UniformOutput',0);
ind = sum([ind{:}],2);

>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8
% Do a quick check to see if we can use a simple mean (equal monitor measurement)
m = y./e.^2;
if all(m==m(1))
    method = 'm';
end

switch lower(method(1))
    case 'm'
        % Simple mean - Remastered to be GPU compatible
        N = accumarray(ind(:),ones(size(x)),[],@sum,NaN);
        xs = accumarray(ind(:),x(:),[],@sum,NaN)./N(:);
        ys = accumarray(ind(:),y(:),[],@sum)./N(:);
        %         es = sqrt(accumarray(ind(:),e(:).^2,[],@sum))./N(:);
        es = accumarray(ind(:),e(:),[],@norm)./N(:);
        if ~isempty(y_fit)
            y_fit_s = accumarray(ind(:),y_fit(:),[],@sum)./N(:);
        else
            y_fit_s = [];
        end
    case 'c'
        % Counts
        % Reconstruct the montor count.
        % For zero error, this is impossible, so we take average monitor value.
        if any(y < 0)
            warning('spec1d:combine:NegativeY','Some Y values are negative. Using method "counts" is not recomended.')
            m = abs(m);
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
        % We have attempted to re-establish sqrt statistics
        es = (accumarray(ind(:),(y(:)./m(:)),[],@sum).^0.5)./ms;
<<<<<<< HEAD
=======
        %         es = accumarray(ind(:),(y(:)./m(:)),[],@norm)./ms;
>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8
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
        es = sqrt(accumarray(ind(:),abs(y(:))./e(:).^2,[],@sum))./ms;%1./(sqrt(accumarray(ind(:),1./e(:).^2,[],@sum)).*accumarray(ind(:),ones(size(x)),[],@sum,NaN));
    otherwise
        error('spec1d:combine:NotValidMethod','%s is not a valid method. See documentation',p.Results.method)
end

% This takes into account the possibility of non-continuous indexing
r = s(1);
r.x = xs(~isnan(xs));
r.y = ys(~isnan(xs));
r.e = es(~isnan(xs));
if ~isempty(y_fit_s)
    r.yfit = y_fit_s(~isnan(xs));
end
s_out = feval(class(r),r);
<<<<<<< HEAD
=======


    function edges = calc_bins(varargin)
        maxx = max(x);
        minx = min(x);
        maximumbins = 62236;
        xrange = maxx - minx;
        hardlimits = false;
        if ~isempty(x) &&  xrange <= 50 && maxx <= flintmax(class(maxx))/2 ...
                && minx >= -flintmax(class(minx))/2
            
            xscale = max(abs(x(:)));
            if nargin == 0
                if xrange > maximumbins
                    % If there'd be more than maximum bins, center them on an appropriate
                    % power of 10 instead.
                    binwidth = 10^ceil(log10(xrange/maximumbins));
                elseif isfloat(x) && eps(xscale) > 1
                    % If a bin width of 1 is effectively zero relative to the magnitude of
                    % the endpoints, use a bigger power of 10.
                    binwidth = 10^ceil(log10(eps(xscale)));
                else
                    % Otherwise bins are centered on integers.
                    binwidth = 1;
                end
            else
                binwidth = varargin{1};
            end
            if ~hardlimits
                minx = binwidth*round(minx/binwidth); % make the edges bin width multiples
                maxx = binwidth*round(maxx/binwidth);
                edges = (floor(minx)-.5*binwidth):binwidth:(ceil(maxx)+.5*binwidth);
            else
                minxi = binwidth*ceil(minx/binwidth)+0.5;
                maxxi = binwidth*floor(maxx/binwidth)-0.5;
                edges = [minx minxi:binwidth:maxxi maxx];
            end
        else
            edges = scottsrule(x,minx,maxx,hardlimits);
        end
    end

    function edges = scottsrule(x, minx, maxx, hardlimits)
        % Scott's normal reference rule
        if ~isfloat(x)
            x = double(x);
        end
        binwidth = 3.5*std(x)/(numel(x)^(1/3));
        if ~hardlimits
            edges = matlab.internal.math.binpicker(minx,maxx,[],binwidth);
        else
            edges = matlab.internal.math.binpickerbl(min(x(:)),max(x(:)),minx,maxx,binwidth);
        end
    end

    function edges = fdrule(x, minx, maxx, hardlimits)
        n = numel(x);
        xrange = max(x(:)) - min(x(:));
        if n > 1
            % Guard against too small an IQR.  This may be because there
            % are some extreme outliers.
            iq = max(datafuniqr(x),double(xrange)/10);
            binwidth = 2 * iq * n^(-1/3);
        else
            binwidth = 1;
        end
        if ~hardlimits
            edges = matlab.internal.math.binpicker(minx,maxx,[],binwidth);
        else
            edges = matlab.internal.math.binpickerbl(min(x(:)),max(x(:)),minx,maxx,binwidth);
        end
    end

    function edges = sturgesrule(x, minx, maxx, hardlimits)
        nbins = max(ceil(log2(numel(x))+1),1);
        if ~hardlimits
            binwidth = (maxx-minx)/nbins;
            if isfinite(binwidth)
                edges = matlab.internal.math.binpicker(minx,maxx,[],binwidth);
            else
                edges = matlab.internal.math.binpicker(minx,maxx,nbins,binwidth);
            end
        else
            edges = linspace(minx,maxx,nbins+1);
        end
    end

    function edges = sqrtrule(x, minx, maxx, hardlimits)
        nbins = max(ceil(sqrt(numel(x))),1);
        if ~hardlimits
            binwidth = (maxx-minx)/nbins;
            if isfinite(binwidth)
                edges = matlab.internal.math.binpicker(minx,maxx,[],binwidth);
            else
                edges = matlab.internal.math.binpicker(minx,maxx,nbins,binwidth);
            end
        else
            edges = linspace(minx,maxx,nbins+1);
        end
    end

    function C = optim_bin(b)
        bb = linspace(b - 0.8*b,b + 0.8*b,10);
        for i = 1:length(bb)
            EDG = calc_bins(bb(i));
            ki = cellfun(@sum,(arrayfun(@(i) (x >= EDG(i)) & (x < EDG(i+1)),1:(length(EDG)-1),'UniformOutput',0)));
            k = mean(ki);
            v = sum((ki - k).^2)/(length(EDG) - 1);
            C(i) = abs((2*k - v)/(bb(i)^2));
        end
        C = interp1(bb,C,bb([1 end]),'pchip');
        C = mean(C);
    end
end

>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8

