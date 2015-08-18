function varargout = sort(s,varargin)
    %% [sout ind ] = sort(s,options)
    % Sort a spec1d object with respect to the x-axis
    % options can be 'ascend' or 'descend'
    % output is sorted spec1d object and the indices if requested.
    %
    % Simon Ward 03/2014, simon.ward@psi.ch
    
    %% Parse inputs
    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired('spec1d',@(x) isa(x,'spec1d'))
    p.addOptional('mode','ascend',@ischar)
    
    p.parse(s,varargin{:})
    
    s_in = p.Results.spec1d;
    mode = p.Results.mode;
    
    %% Do the sort
    for i = 1:length(s_in)
        s_out(i) = s_in(i);
        [s_out(i).x, ind] = sort(s_out(i).x,mode);
        s_out(i).y = s_out(i).y(ind);
        s_out(i).e = s_out(i).e(ind);
        if ~isempty(s_out(i).yfit)
            s_out(i).yfit = s_out(i).yfit(ind);
        end
    end
    
    %% Create outputs
    varargout{1} = s_out;
    
    if nargout > 1
        varargout{2} = ind;
    end