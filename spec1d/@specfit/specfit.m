classdef specfit
    %SPECFIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        pvals
    end
    
    properties
        evals
        func = '';
        pnames = {};
        userdata
        notfixed = true(0,0);
        parser
    end
    
    properties (SetAccess=protected)
        fitdata = mergestruct(optimset(),struct('CriteriaFcn',[],'JacobFcn',[],'Bounds',[]));
    end
    
    properties (Dependent)
        chisq
        rsq
    end
    
    properties(Hidden=true)
        specID
        p_vals
    end
    
    properties (SetAccess=immutable,Hidden=true)
        ident= java.rmi.server.UID();
    end
    
    methods
        function obj = specfit(varargin)
            if nargin == 0
                return
            elseif nargin == 3
                obj.func = varargin{1};
                obj.p_vals = varargin{2};
                obj.evals = nan(size(obj.pvals));
                obj.notfixed = varargin{3};
            elseif nargin >= 4
                % Fill out the basics for creating a fit object.
                obj.p_vals = varargin{1};
                obj.evals = varargin{2};
                obj.func = varargin{3};
                obj.notfixed = varargin{4};
                if nargin > 4
                    % If asigning the results we asign a spec1d object.
                    if isa(varargin{5},'spec1d')
                        obj.specID = varargin{5}.ident;
                    elseif isa(varargin{5},'java.rmi.server.UID')
                        obj.specID = varargin{5};
                    else
                        obj.pnames = varargin{5};
                    end
                    % Add pnames if required.
                    if nargin > 5
                        obj.pnames = varargin{6};
                    end
                end
                % Structural asignment just because.
            elseif nargin == 1
                if isstruct(varargin{1})
                    f = fieldnames(varargin{1});
                    for i = 1:length(f)
                        try
                            obj.(f{i}) = varargin{1}.(f{i});
                        end
                    end
                elseif isa(varargin{1},'specfit')
                    obj = varargin{1}.copy;
                end
            end
        end
        
        function y = evaluate(obj,x,varargin)
            if nargin == 1
                x = findspec1d(obj);
                if isempty(x)
                    y = [];
                    return
                end
            end
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
        
        function value = get.chisq(obj)
            value = 0;
            if isempty(obj.pvals) || isempty(obj.notfixed)
                return
            end
            s = findspec1d(obj);
            if isempty(s)
                return
            end
            v = length(s.y)-sum(logical(obj.notfixed));
            value = sum(((s.y-s.yfit)./s.e).^2 )/v;
        end
        
        function p = get.pvals(obj)
            if ~isstruct(obj.p_vals)
                p = obj.p_vals;
            else
                if ~isempty(obj.parser)
                    p = [];
                    for i = 1:length(obj.parser) 
                        p = vertcat(p,reshape(eval(sprintf('obj.p_vals.%s',obj.parser{i})),[],1));
                    end
                else
                    p = 'This is a structure, please define a parser';
                end
            end
        end
        
        function obj = set.pvals(obj,vals)
            if ~isempty(obj.parser)
                offset = 1;
                for i = 1:length(obj.parser)
                    l = length(eval(sprintf('obj.p_vals.%s',obj.parser{i})));
                    eval(sprintf('obj.p_vals.%s = vals(%i:%i)',obj.parser{i},offset,offset+l-1));
                    offset = offset + l;
                end
            else
                obj.p_vals = vals;
            end
        end
        
        function value = get.rsq(obj)
            value = NaN;
            if isempty(obj.pvals) || isempty(obj.notfixed)
                return
            end
            s = findspec1d(obj);
            if isempty(s)
                return
            end
            r = corrcoef([s.y(:),s.yfit(:)]);
            value = r(1,2).^2;
            % This is the reduced r2....
            %             value = 1-(1-r2)*((length(s.y)-1)/(length(s.y)-sum(obj.notfixed ~= 0)-1));
        end
        
        function st = toStruct(obj)
            st = struct();
            f = properties(obj(1));
            for i = 1:length(f)
                st.(f{i}) = obj.(f{i});
            end
        end
        
        function obj = setFitdata(obj,varargin)
            f = fieldnames(obj.fitdata);
            for i = 1:2:length(varargin)
               if any(strcmp(varargin{i},f))
                   try
                       obj.fitdata.(varargin{i}) = varargin{i+1};
                   catch
                      warning('Can not set property: %s',varargin{i}) 
                   end
               end
            end
        end
        
        function len = lengthP(obj)
            if isempty(obj.parser)
                len = length(obj.pvals);
            else
                len = 0;
                for i = 1:length(obj.parser)
                    l = length(eval(sprintf('obj.p_vals.%s',obj.parser{i})));
                    len = len + l;
                end
            end
        end
    end
    
    methods (Hidden=true)
        function obj_r = copy(obj,obj2)
            obj_r = feval(class(obj));
            if nargin==2
                obj = obj2;
            end
            f = properties(obj);
            for i = 1:length(f)
                try
                    obj_r.(f{i}) = obj.(f{i});
                end
            end
            obj_r.p_vals = obj.p_vals;
            obj_r.specID = obj.specID;
        end
    end
end

