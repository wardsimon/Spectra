function out =  toll(fun,varargin)
    % Check for different types of function.
    % is a string, make a function
    if isa(fun,'char')
        fun = str2func(fun);
        % is a anonymous function, continue
    elseif ~isa(fun,'function_handle')
        % is a function, built in or in the path?
        if ~any(exist(fun)==[2 5])
            error('unknown function')
        end
        % is not a function, do numerics
    elseif isnumeric(fun)
        val = varargin{1};
        tol = varargin{2};
        out = all(fun(:)>=val(:)-tol & fun(:)<=val(:)+tol);
        return
    end
    
    if length(varargin)==3
        % Do a simple evaluation
        if ~iscell(varargin{1})
            p = {varargin{1}};
        else
            p = varargin{1};
        end
        out = feval(fun,p{:});
        val = varargin{2};
        tol = varargin{3};
        out = all(out(:)>=val(:)-tol & out(:)<=val(:)+tol);
    else
        test='';
        for i=1:length(varargin)
            if isa(varargin{i},'char');
                test = varargin{i};
            end
        end
        if isempty(test)
            error('No test given')
        end
        switch test
            case 'abspm'
                temp = 1:length(varargin);
                temp(i)=[];
                temp2 = temp((end-1):(end));
                temp((end-1):(end))=[];
                % We are given a crazy test..
                out = feval(fun,varargin{temp});
                val = varargin{temp2(1)};
                tol = varargin{temp2(2)};
                % Need to write tests....
                out = all(out>=val-tol & out<=val+tol);
        end
    end
end