function obj = setfitdata(obj,fd,varargin)
%SETFITDATA Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(obj)
    if length(obj)==length(fd)
        if isa(fd(i),'specfit')
            fd(i).specID = obj(i).ident;
            obj(i).fitdata = fd(i);
        elseif isstruct(fd(i))
            temp = specfit(fd(i));
            temp.specID = obj(i).ident;
            obj(i).fitdata = temp;
        end
    else
        if isa(fd,'specfit')
            fd.specID = obj(i).ident;
            obj(i).fitdata = fd;
        elseif isstruct(fd)
            temp = specfit(fd);
            temp.specID = obj(i).ident;
            obj(i).fitdata = temp;
        else
            temp = specfit(fd,varargin{:});
            temp.specID = obj(i).ident;
            obj(i).fitdata = temp;
        end
    end
end
end

