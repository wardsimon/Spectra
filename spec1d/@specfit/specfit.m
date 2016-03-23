classdef specfit < handle
    %SPECFIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetObservable)
        pvals
    end
    
    properties
        evals
        func
        pnames
        chisq
        rsq
        notfixed
    end

    properties(Hidden=true)
        specID
    end
    
    properties (SetAccess=protected,Hidden=true)
        listners
    end
    
    properties (SetAccess=immutable,Hidden=true)
        ident;
    end
    
    methods
        function obj = specfit(varargin)
            obj.ident = java.rmi.server.UID();
            if nargin == 0
                return
            elseif nargin == 6
                obj.pvals = varargin{1};
                obj.evals = varargin{2};
                obj.func = varargin{3};
                obj.pnames = varargin{4};
                
            elseif nargin == 1
                if isstruct(varargin{1})
                    f = fieldnames(varargin{1});
                    for i = 1:length(f)
                        try
                            obj.(f{i}) = varargin{1}.(f{i});
                        end
                    end
                end
            end
            %             obj.listners.pvals = addlistener(obj,'pvals','PostSet',@obj.handlePropertyEvents);
        end
        
        function y = evaluate(obj,x,varargin)
            if isa(x,'spec1d')
                for i = 1:length(x)
                    if isempty(varargin)
                        y(i) = feval(obj.func,x(i).x,obj.pvals);
                    else
                        y(i) = feval(obj.func,x(i).x,obj.pvals,varargin{:});
                    end
                end
            else
                if isempty(varargin)
                    y = feval(obj.func,x,obj.pvals);
                else
                    y = feval(obj.func,x,obj.pvals,varargin{:});
                end
            end
        end 
    end
    
    methods (Hidden=true)
        function obj2 = equals(obj)
            for i = 1:length(obj);
                obj2(i) = obj(i).copy;
            end
        end
        
        function obj2 = copy(obj)
            obj2 = specfit();
            f = properties(obj);
            for i = 1:length(f)
                try
                    obj2.(f{i}) = obj.(f{i});
                end
            end
        end
    end
end

