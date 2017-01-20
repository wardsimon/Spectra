classdef spec1d
        
    properties (SetAccess=protected)
        x
        y
        e
        yfit
        fitdata
    end
    
    properties
        userdata
    end
    
<<<<<<< HEAD
    if length(a) > 1 && isstruct(a)
        si = size(a);
        a = a(:);
        for i = 1:length(a)
            [a(i).x, ind] = sort(a(i).x);
            a(i).y = a(i).y(ind);
            a(i).e = a(i).e(ind);
            if ~isempty(a(i).yfit)
                a(i).yfit = a(i).yfit(ind);
            end
            s_in(i) = spec1d(a(i));
        end
        s = reshape(s_in,si);
        return
    else
        p.parse(a,varargin{:});
        s_in = p.Results;
=======
    properties (Hidden=true)
        datafile
        x_label
        y_label
>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8
    end
    
    properties (SetAccess=immutable,Hidden=true)
        ident;
    end
        
    methods
        function s = spec1d(varargin)
            s.ident = java.rmi.server.UID();
            if nargin == 0
                s_in = constructor();
                f = fieldnames(s_in);
                for i = 1:length(f)
                    s.(f{i}) = s_in.(f{i});
                end
                s.fitdata.specID = s.ident;
                return
            else
%                 s = spec1d;
                s_in = constructor(varargin{:});
                f = fieldnames(s_in);
                for i = 1:length(f)
                    s.(f{i}) = s_in.(f{i});
                end
                validate(s);
                s.fitdata.specID = s.ident;
            end
        end
                
        function disp(s)
            if length(s) == 1
                if isempty(s.datafile)
                    fprintf('Spec1d object with %i datapoint/s.\n',numel(s.x))
                else
                    fprintf('Spec1d object: %s with %i datapoint/s.\n',s.datafile,numel(s.x))
                end
            else
                fprintf('Spec1d object of size [%i %i]\n',size(s,1),size(s,2))
            end
        end
        
        function s_out = pushToGraphicsCard(s)
            if ~checkGpuMemory(s)
                warning('spec1d:pushToGraphicsCard','There is not enough space on the graphics card or it''s not functioning. Aborting...')
                s_out = s;
                return
            end
            s_out = spec1d;
            p = properties(s);
            for i = 1:length(p)
                if any(strcmp(p{i},{'x','y','e'}))
                    s_out.(p{i}) = gpuArray(s.(p{i}));
                else
                    s_out.(p{i}) = s.(p{i});
                end
            end
        end
        
        function obj = equals(obj)
            obj = obj.copy; % We do this so that the identity is always different.
        end
        
    end
    
    methods (Hidden=true)
       out = isempty(obj);
       out = findID(obj,ID);
    end
    
    methods (Access=protected)
        names = fieldnames_s(obj); % Super/sub classes can use this.
    end
end