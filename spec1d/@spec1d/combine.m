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
% 'relative'    : Histogram bining of size 'toll', averages are taken in
%                 this range. The minimum gap before averaging is greater
%                 than 'toll'
% 'absolute'	: The starting point of the histogram is continuously shifted
%                 to the data points. The absolute bin between spectra points
%                 is always greater then 'toll'.
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

p = inputParser;
p.CaseSensitive = false;
p.addRequired('toll',@(x) isnumeric(x) && isreal(x));
p.addRequired('s_in',@(x) isa(x,'spec1d'));
p.addOptional('add_s',[],@(x) isa(x,'spec1d'))
p.addOptional('indexing','relative',@ischar)
p.addParameter('method','counts',@ischar);
p.parse(toll,varargin{:});
if ~isempty(p.Results.add_s(:))
    s = [p.Results.s_in(:) p.Results.add_s(:)];
else
    s = p.Results.s_in(:);
end
bin = p.Results.toll;
method = p.Results.method;
indexing = p.Results.indexing;

x = [s.x];
[x, ind ] = sort(x(:));
y = [s.y]; y = y(ind);
e = [s.e]; e = e(ind);

gpuCompute = 0;
if isa(x,'gpuArray')
    gpuCompute = 1;
    if ~strcmp(method,{'counts','weight'})
        warning('spec1d:combine:GpuUnsuportedMethod','%s is not a valid GPU method. Switching to counts',p.Results.method)
        method = 'counts';
    end
end

switch lower(indexing(1))
    case 'r'
        % Relative
        ind = [1; ceil(cumsum(diff(x(:)))/bin)]; % This is a faster way of [~,~,ind] = histcounts(x(:),min(x):bin:max(x));
        % - Using accumarray(~,~,[],@sum,NaN) and filtering NaNs we can remove these time consuming statements -
        %         continuous = diff(ind);
        %         if any(continuous > 1) % This makes the vector monotonically increasing.
        %             while any(continuous > 1)
        %                 start_index = find(continuous>1,1,'first')+1;
        %                 ind(start_index:end) =  ind(start_index-1) + ind(start_index:end) - ind(start_index) + 1;
        %                 continuous = diff(ind);
        %             end
        %         end
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
            if sum(diff(xtemp)) > bin
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
        % Simple mean - This is not currently supported by the GPU due to
        % @mean - As of R2016a (pre-release)
        xs = accumarray(ind(:),x(:),[],@mean,NaN);
        ys = accumarray(ind(:),y(:),[],@mean);
        es = accumarray(ind(:),e(:),[],@(x) sqrt(sum(x.^2))/length(x));
    case 'c'
        % Counts
        % Reconstruct the montor count.
        % For zero error, this is impossible, so we take average monitor value.
        if gpuCompute
            m = zeros(size(e),'gpuArray');
        else
            m = zeros(size(e));
        end
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
        es = (accumarray(ind(:),(e(:).*m(:)).^2,[],@sum).^0.5)./ms;
    case 'w'
        % Weight
        ms = accumarray(ind(:),1./e(:).^2,[],@sum);
        ms(ms < eps) = eps; % Check for possibility of rounding errors
        xs = accumarray(ind(:),x(:)./e(:).^2,[],@sum,NaN)./ms;
        ys = accumarray(ind(:),y(:)./e(:).^2,[],@sum)./ms;
        es = 1./sqrt(accumarray(ind(:),1./e(:).^2,[],@sum));
    otherwise
        error('spec1d:combine:NotValidMethod','%s is not a valid method. See documentation',p.Results.method)
end

s_out = spec1d(xs(~isnan(xs)),ys(~isnan(xs)),es(~isnan(xs))); % This takes into account the possibility of non-continuous indexing

