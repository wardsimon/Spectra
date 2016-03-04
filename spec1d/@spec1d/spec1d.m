classdef spec1d
    properties
        x
        y
        e
        yfit
        datafile
        x_label
        y_label
        userdata
    end
    methods
        function s = spec1d(varargin)
            if nargin == 0
                return
            else
                s = spec1d;
                p = properties(s);
                s_in = constructor(varargin{:});
                for i = 1:length(p)
                    s.(p{i}) = s_in.(p{i});
                end
                validate(s);
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
    end
end