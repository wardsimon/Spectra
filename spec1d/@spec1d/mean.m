function varargout = mean(varargin)
%
% function [mean error] = mean(s1..sn,'method','mean')
%
% @SPEC1D/mean function to give the mean and optionally error of each spectrum s1...sn.
% Optional arguments are 'method' and a choice of:
% 'mean' - A simple mean (default)
% 'counts' - Weighted by counts and error
% 'weight' - Weighted by error
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%
s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin(s_ind) = [];

p = inputParser;
p.CaseSensitive = false;
p.addRequired('s_in',@iscell);
p.addParameter('method','mean',@ischar);
p.parse(s,varargin{:});

method = p.Results.method;
s = [p.Results.s_in{:}];

m_out = zeros(size(s));
e_out = m_out;

for i = 1:length(s)
    y = s(i).y(:);
    e = s(i).e(:);
    
    gpuCompute = 0;
    if isa(y,'gpuArray')
        gpuCompute = 1;
        if ~strcmp(method,'counts')
            warning('spec1d:mean:GpuUnsuportedMethod','%s is not a valid method. Switching to counts',p.Results.method)
            method = 'counts';
        end
    end
    
    ind = ones(size(y));
    switch method
        case 'mean'
            % Simple mean
            ys = accumarray(ind(:),y(:),[],@mean);
            es = accumarray(ind(:),e(:),[],@(x) sqrt(sum(x.^2))/length(x));
        case 'counts'
            % Counts
            if gpuCompute
                m = zeros(size(e),'gpuArray');
            else
                m = zeros(size(e));
            end
            m(e~=0) = y(e~=0)./e(e~=0).^2;
            if any(e==0)
                m(e==0) = mean(m(e~=0));
            end
            ms = accumarray(ind(:),m(:),[],@sum);
            ys = accumarray(ind(:),y(:).*m(:),[],@sum)./ms;
            es = (accumarray(ind(:),(e(:).*m(:)).^2,[],@sum).^0.5)./ms;
        case 'weight'
            % Weight
            ms = accumarray(ind(:),1./e(:),[],@(x) max(sum(x),eps));
            ys = accumarray(ind(:),y(:)./e(:),[],@sum)./ms;
            es = accumarray(ind(:),(1./e(:)).^2,[],@(x) sqrt(sum(x)))./ms;
        otherwise
            error('spec1d:mean:NotValidMethod','%s is not a valid method. See documentation',p.Results.method)
    end
    m_out(i) = ys;
    e_out(i) = es;
end

varargout{1} = m_out;
if nargout == 2
    varargout{2} = e_out;
end
