function s = findspec1d( obj )
%FINDSPEC1D Summary of this function goes here
%   Detailed explanation goes here
v = evalin('base','who');
id = obj.specID;
s = [];
for i = 1:length(v)
    temp = evalin('base',sprintf('isa(%s,''spec1d'')',v{i})); 
    if temp
        s = evalin('base',sprintf('%s',v{i})); 
        if any(s.findID(id))
            return
        end
    end; 
end

end

