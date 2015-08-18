function varargout = get(s,varargin)
    
    if isa(s,'spec1d')
        f_names = fieldnames(s(1));
    else
        error('No spec1d object specified')
    end
    
    for i = 1:length(varargin)
        if ~any(strcmpi(varargin{i},f_names))
            error('%s is not a valid field',varargin{i})
        end
    end
    for i = 1:length(s)
        varargout{i} = zeros(length(s(i).x),length(varargin));
        for j = 1:length(varargin)
            if nargin-1 == nargout && length(s)==1 
                varargout{j} = s(i).(varargin{j});
            else
                varargout{i}(:,j) = s(i).(varargin{j});
            end
        end
    end