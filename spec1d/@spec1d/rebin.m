function s_out = rebin(s,dx,varargin)
% function s=rebin(s, dx)
%
% SPED1D/REBIN Rebins a spec1d spectrum.
%
% dx = scalar: specifies bin widths.
% dx = vector: specifies the bin centers.
% dx = sped1d: uses x-values of dx as bin centers
% dx = string: Use a binning algorithm described in the documentation histcounts.
%
% The width of boundary bins is twice the distance to neighbour bin.
%
% Optional parameter pairs.
% method 'average' : New x values are average within bins
%                    (weighted by error). 'Default'
%        'force'   : New x values are as specified in dx
%        'interp'  : New x-values as specified by dx, y and e
%                    interpolated from 'average'-x,
% byError 0        : Data is binned by x. 'Default'
%         1        : Data is binned to a constant error. Only 'method',
%                    average is supported.
% Simon Ward 28/01/2016 - simon.ward@psi.ch

p = inputParser;
p.addRequired('s',@(x)isa(x,'spec1d'));
p.addRequired('dx',@(x) (isnumeric(x) && isreal(x)) || (isa(x,'spec1d') && length(x)==1) || ischar(x))
p.addParameter('method','average',@ischar)
p.addParameter('byError',0,@isnumeric)

p.parse(s,dx,varargin{:})

s = p.Results.s;
dx = p.Results.dx;
method = p.Results.method;
if ischar(dx) && verLessThan('matlab','8.4')
    error('spec1d:rebin:MatlabVerionError','This option needs the function histcounts which is not available on your MATLAB version.')
end

for i = 1:length(s)
    x = get(s(i),'x');
    y = get(s(i),'y');
    e = get(s(i),'e');
    yfit = get(s(i),'yfit');
    
    if p.Results.byError
        if any(strcmpi(method(1),{'f','i'}))
            warning('spec1d:rebin:byErrorMethodNotSupported','%s is not a method when using byError. Switching to "average".',method)
            method = 'a';
        end
        [N,ind] = error_index(e,dx);
    else
        % Deal with dx options
        if isa(dx,'spec1d');
            % This is the spec1d case
            centers = get(dx,'x');
            d = diff(centers)/2;
            edges = [centers(1)-d(1); centers(1:end-1)+d; centers(end)+d(end)];
            if verLessThan('matlab','8.4')
                [N,ind] = histc(x(:),edges);
            else
                [N,~,ind] = histcounts(x(:),edges);
            end
        elseif ischar(dx)
            % This is the histcounts binning algorithm case
            try
                [N,edges,ind] = histcounts(x(:),'BinMethod',dx); % Matlab supplied algorithm
            catch
                error('spec1d:rebin:NotValidAlgorithm','%s is not a valid binning algorithm for histcounts. See documentation.',dx)
            end
        else
            % This is the numeric case
            if length(dx) == 1;
                % This is the scalar case
                if verLessThan('matlab','8.4')
                    edges = linspace(min(x),max(x),round((max(x)-min(x))/dx));
                    [N, ind] = histc(x(:),edges);
                else
                    [N,edges,ind] = histcounts(x(:),'BinWidth',dx);
                end
            else
                % This is the vecor case
                centers = dx;
                d = diff(centers)/2;
                edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];
                if verLessThan('matlab','8.4')
                    [N, ind] = histc(x(:),edges);
                else
                    [N,~,ind] = histcounts(x(:),edges);
                end
            end
        end
    end
    % Take care if edges does not cover data range. Remove the uncovered
    % data.
    if any(ind == 0)
        warning('spec1d:rebin:InvalidBinRange','Supplied bin range does not cover data. Removing points.')
        x = x(ind ~= 0);
        y = y(ind ~= 0);
        e = e(ind ~= 0);
        if ~isempty(yfit)
            yfit = yfit(ind ~= 0);
        end
        ind = ind(ind ~= 0);
    end
    
    % Perform the binning
    switch lower(method(1))
        case 'a'
            % Average
            ms = accumarray(ind(:),1./e(:).^2,[],@sum);
            ms(ms < eps) = eps; % Check for possibility of rounding errors
            xs = accumarray(ind(:),x(:)./e(:).^2,[],@sum,NaN)./ms;
            ys = accumarray(ind(:),y(:)./e(:).^2,[],@sum)./ms;
            es = sqrt(accumarray(ind(:),e(:).^2,[],@sum))./N(:);
            if ~isempty(yfit)
                yfits = accumarray(ind(:),yfit(:)./e(:).^2,[],@sum)./ms;
            end
        case 'i'
            % Interpolate
            xs = edges(1:(end-1)) + diff(edges)/2;
            temp = interpolate(s(i),xs);
            xs = get(temp,'x');
            ys = get(temp,'y');
            es = get(temp,'e');
            if ~isempty(yfit)
                yfits = get(temp,'yfit');
            end
        case 'f'
            % Force xs, convert bin edges to centers
            xs = edges(1:(end-1)) + diff(edges)/2;
            ys = accumarray(ind(:),y(:),[],@sum)./N(:);
            es = sqrt(accumarray(ind(:),e(:).^2,[],@sum))./N(:);
            if ~isempty(yfit)
                yfits = accumarray(ind(:),yfit(:),[],@sum)./N(:);
            end
        otherwise
            error('spec1d:rebin:NotValidMethod','%s is not a valid method. See documentation',p.Results.method)
    end
    
    % Make output object
    r = s(i);
    r.x = xs(~isnan(xs));
    r.y = ys(~isnan(xs));
    r.e = es(~isnan(xs));
    if ~isempty(r.yfit)
        r.yfit = yfits(~isnan(xs));
    end
    s_out = spec1d(r);
    
end

end

function [N, ind] = error_index(e,bin)
e2 = e(:).^2;

N = nan(size(e2));
ind = zeros(size(e2));

ind(1) = 1;
xtemp = e2(1);
j = 1;
for k = 2:length(e)
    if k == 2
        if sqrt(sum(xtemp))/length(xtemp) < bin
            ind(k) = ind(k-1)+1;
            N(j) = length(xtemp);
            j = j+1;
            xtemp = [];
        else
            ind(k) = ind(k-1);
        end
    else
        if isempty(xtemp)
            xtemp = e2(k);
        else
            xtemp = [xtemp e2(k)]; %#ok<AGROW>
        end
        if sqrt(sum(xtemp))/length(xtemp) < bin
            ind(k) = ind(k-1)+1;
            N(j) = length(xtemp);
            j = j+1;
            xtemp = [];
        else
            ind(k) = ind(k-1);
        end
    end
end
N(isnan(N)) = [];

if sum(max(ind)) == length(N)+1
    N(end+1) = 1;
end
end