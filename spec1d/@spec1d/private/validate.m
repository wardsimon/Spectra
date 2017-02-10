function  validate( s_in )
% Validates a spec1d object to make sure all fields are correct.
%
% VALIDATE(s_in)
%
% Input:
%
% s_in  A spec1d object which needs to be validated.
%

p  = inputParser;
p.addParamValue('x',@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'x'))
p.addParamValue('y',@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'y'))
p.addParamValue('e',@(x) validateattributes(x,{'numeric','gpuArray'},{'nonnegative','vector','real','finite','nonnan'},mfilename,'e'))
p.addParamValue('yfit',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'real','finite','nonnan','2d'},mfilename,'yfit'))
p.addParamValue('fitdata',specfit(),@(x)isa(x,'specfit'))
p.addParamValue('x_label','',@(x)validateattributes(x,{'char'},{'2d'},mfilename,'x_label'))
p.addParamValue('y_label','',@(x)validateattributes(x,{'char'},{'2d'},mfilename,'y_label'))
p.addParamValue('datafile','',@(x)validateattributes(x,{'char'},{'2d'},mfilename,'datafile'))
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

