function sout = set(s,varargin)
    
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = false;
    p.addRequired('spec1d',@(x) isa(x,'spec1d'));
    names = fieldnames(spec1d());
    for i = 1:length(names)
        p.addParamValue(names{i},NaN); %#ok<NVREPL>
    end
    p.parse(s,varargin{:});
    
    s = p.Results.spec1d;
    
    opt = p.Results;
    opt = rmfield(opt,'spec1d');
   
    sout = s.copy;
    for i = 1:length(names)
        if ~isnan(opt.(names{i}))
            sout.(names{i}) = opt.(names{i});
        end
    end
    sout = spec1d(sout);