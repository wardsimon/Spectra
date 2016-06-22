classdef specfit
    %SPECFIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pvals
        evals
        func = '';
        pnames = {};
        userdata
        notfixed = true(0,0);
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
                obj.pvals = varargin{2};
                obj.evals = nan(size(obj.pvals));
                obj.notfixed = varargin{3};
            elseif nargin >= 4
                % Fill out the basics for creating a fit object.
                obj.pvals = varargin{1};
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
               if strcmp(varargin{i},f)
                   try
                       obj.(varargin{i}) = varargin{i+1};
                   catch
                      warning('Can not set property: %s',varargin{i}) 
                   end
               end
            end
        end
        
    end
    
    methods (Hidden=true)
        
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

