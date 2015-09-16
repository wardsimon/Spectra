function s_out = power(varargin)
    p = inputParser;
    p.addRequired('spec1d',@(x) isa(x,'spec1d'))
    p.addRequired('power',@isnumeric)
    
    p.parse(varargin{:})
    
    results = p.Results;
    
    s = results.spec1d;
    p = results.power;
    
    for i = 1:length(s)
        s_out(i) = s(i);
        s_out(i).y = s_out(i).y.^p;
        s_out(i).e = p*s_out(i).y .* (s_out(i).e./s(i).y);
    end
    
        