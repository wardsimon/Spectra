function s = set(s,varargin)
    
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = false;
    p.addRequired('spec1d',@(x) isa(x,'spec1d'));
    p.addParamValue('x',[],@(x) isnumeric(x) & isreal(x))
    p.addParamValue('y',[],@(x) isnumeric(x) & isreal(x))
    p.addParamValue('e',[],@(x) isnumeric(x) & isreal(x))
    p.addParamValue('yfit',[],@(x) isnumeric(x) & isreal(x))
    p.addParamValue('x_label',[],@(x) ischar(x) | isempty(x))
    p.addParamValue('y_label',[],@(x) ischar(x) | isempty(x))
    p.addParamValue('datafile',[],@(x) ischar(x) | isempty(x))
    p.parse(s,varargin{:});
    
    s = p.Results.spec1d;
    
    opt = p.Results;
    opt.x = opt.x(:);
    opt.y = opt.y(:);
    opt.e = opt.e(:);
    opt.yfit = opt.yfit(:);
    
    if isempty(s)
        error('No spec1d object specified')
    end
    
    for i = 1:length(s)
        for j = 1:length(p.Parameters)
            if ~any(strcmp(p.Parameters{j},p.UsingDefaults)) && ~strcmp(p.Parameters{j},'spec1d')
                s(i).(p.Parameters{j}) = opt.(p.Parameters{j});
            end
        end
        if any(diff([length(s(i).x) length(s(i).y) length(s(i).e)]))
            warning('x, y and e are not of the same length!\nLength x:\t%i\nLength y:\t%i\nLength e:\t%i',length(s(i).x),length(s(i).y),length(s(i).e))
        end
    end
    
