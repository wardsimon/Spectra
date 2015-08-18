function varargout = inde(varargin)
    
    narginchk(2,inf)
    ind=varargin{end};
    nargoutchk(1, nargin-1)
    
    for i=1:(nargin-1)
        try
            varargout{i}=varargin{i}(ind);
        catch err
            warning(sprintf('inde can not complete on input %i:\n%s',i,err.message))
            varargout{i}=NaN;
        end
    end
end