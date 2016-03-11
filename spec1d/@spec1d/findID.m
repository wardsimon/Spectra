function [ output_args ] = findID(obj,ID)
%FINDID Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(obj)
    output_args(i) = obj(i).ident.equals(ID);
end

end

