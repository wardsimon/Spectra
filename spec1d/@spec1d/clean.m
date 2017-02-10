function obj = clean(obj)
%CLEAN Summary of this function goes here
%   Detailed explanation goes here

        % Clean up e
        ind = ((obj.e == 0 | isinf(obj.e) | isnan(obj.e)) | ...
            (isinf(obj.y) | isnan(obj.y)) | ...
            (isinf(obj.x) | isnan(obj.x)));
        
        obj.x(ind) = [];
        obj.y(ind) = [];
        obj.e(ind) = [];
        
        if ~isempty(obj.yfit)
            obj.yfit(ind) = [];
        end
        
        if isfield(obj.userdata,'rind')
            if ~isempty(obj.userdata.rind)
                obj.userdata.rind(ind) = [];
            end
        end

        obj = spec1d(obj);
end

