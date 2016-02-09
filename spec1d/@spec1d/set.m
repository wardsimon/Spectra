function s = set(s,varargin)
    
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = false;
    p.addRequired('spec1d',@(x) isa(x,'spec1d'));
    names = fieldnames(spec1d());
    for i = 1:length(names)
        p.addParamValue(names{i},[]); %#ok<NVREPL>
    end
    p.parse(s,varargin{:});
    
    s = p.Results.spec1d;
    
    opt = p.Results;
    opt = rmfield(opt,'spec1d');
   
    s = setstructfields(s,opt);
    s = spec1d(s);