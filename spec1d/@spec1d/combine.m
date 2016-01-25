function s_out = combine(toll,varargin)
%
% function r=combine(method,toll,s1,s2,....sn)
%
% @SPEC1D/COMBINE function to combine two or more spectra. 
%
% If the x values of two points differ by less
% than tolerance toll, then the points are combined.
%
% Depending on method, points are combined as
% 'mean'		: Simple means for x and y, errors are averaged in quadrature.
% 'counts'	: Restablishes normalisation and original counts assuming 
%				  square-root statistics, is correct for normalised counts
% 'weight'	: Weights to inverse error. For more general data.
% Default is 'counts'
%
% s1,s2,... can be single spectra or arrays of spectra.
%
% Example: 
% Combine s1,s2 and s3 when x values differ by less than 0.01.
% >r=combine('mean',0.01,s1,s2,s3)
%
% DFM 1.4.98, HMR, NBC, BHL 20.11.2000

p = inputParser;
p.addRequired('toll',@(x) isnumeric(x) && isreal(x));
p.addRequired('s_in',@(x) isa(x,'spec1d'));
p.addOptional('add_s',[],@(x) isa(x,'spec1d'))
p.addParameter('method','counts',@ischar);
p.parse(toll,varargin{:});
s = [p.Results.s_in(:) p.Results.add_s(:)];
bin = p.Results.toll;

x = [s.x]; [x, ind ] = sort(x(:));
y = [s.y]; y = y(ind);
e = [s.e]; e = e(ind);


ind = [1; ceil(cumsum(diff(x(:)))/bin)];

% end

switch p.Results.method
    case 'mean'
        % Simple mean
        xs = accumarray(ind(:),x(:),[],@mean);
        ys = accumarray(ind(:),y(:),[],@mean);
        es = accumarray(ind(:),e(:),[],@(x) sqrt(sum(x.^2))/length(x));
    case 'counts'
        % Counts
        m = zeros(size(e));
        m(e~=0) = y(e~=0)./e(e~=0).^2;
        if any(e==0)
            m(e==0) = mean(m(e~=0));
        end
        ms = accumarray(ind(:),m(:),[],@(x) max(sum(x),eps));
        xs = accumarray(ind(:),x(:).*m(:),[],@sum)./ms;
        ys = accumarray(ind(:),y(:).*m(:),[],@sum)./ms;
        es = accumarray(ind(:),(e(:).*m(:)).^2,[],@(x) sqrt(sum(x)))./ms;
    case 'weight'
        % Weight
        ms = accumarray(ind(:),1./e(:),[],@(x) max(sum(x),eps));
        xs = accumarray(ind(:),x(:)./e(:),[],@sum)./ms;
        ys = accumarray(ind(:),y(:)./e(:),[],@sum)./ms;
        es = accumarray(ind(:),(1./e(:)).^2,[],@(x) sqrt(sum(x)))./ms;
    otherwise
        error('spec1d:combine:NotValidMethod','%s is not a valid method. See documentation',p.Results.method)
end

s_out = s(1);
s_out.x = xs;
s_out.y = ys;
s_out.e = es;

