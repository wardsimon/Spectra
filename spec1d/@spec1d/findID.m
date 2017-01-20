function [ output_args ] = findID(obj,ID)
%FINDID Summary of this function goes here
%   Detailed explanation goes here

k = 1;
for i = 1:length(obj)
    for j = 1:length(ID)
        output_args(k) = obj(i).ident.equals(ID(j));
        k = k+1;
    end
end

end

