function  validate( s_in )
%VALIDATE Summary of this function goes here
%   Detailed explanation goes here
p  = inputParser;
p.addParamValue('x',@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'x'))
p.addParamValue('y',@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'y'))
p.addParamValue('e',@(x) validateattributes(x,{'numeric','gpuArray'},{'nonnegative','vector','real','finite','nonnan'},mfilename,'e'))
p.addParamValue('yfit',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'real','finite','nonnan','2d'},mfilename,'yfit'))
p.addParamValue('x_label',[],@(x) isempty(x) ||validateattributes(x,{'char'},{},mfilename,'x_label'))
p.addParamValue('y_label',[],@(x) isempty(x) ||validateattributes(x,{'char'},{},mfilename,'y_label'))
p.addParamValue('datafile',[],@(x) isempty(x) ||validateattributes(x,{'char'},{},mfilename,'datafile'))
p.addParamValue('userdata',[])

a = struct;
pr = properties(s_in);
for i = 1:length(pr)
    a.(pr{i}) = s_in.(pr{i});
end

p.parse(a);
s_in_p = p.Results;

    if any(diff([length(s_in_p.x) length(s_in_p.y) length(s_in_p.e)]))
        warning('x, y and e are not of the same length!\nLength x:\t%i\nLength y:\t%i\nLength e:\t%i',length(s_in_p.x),length(s_in_p.y),length(s_in_p.e))
    end

end

