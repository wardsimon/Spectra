function s = findspec1d( obj )
%FINDSPEC1D Summary of this function goes here
%   Detailed explanation goes here
v = evalin('base','who');
id = obj.specID;
s = [];
for i = 1:length(v)
    temp = evalin('base',sprintf('isa(%s,''spec1d'')',v{i})); 
    if temp
        s_in = evalin('base',sprintf('%s',v{i})); 
        for j = 1:length(s_in)
            if s_in(j).findID(id)
                s = s_in(j);
                return
            end
        end
    end; 
end

end

