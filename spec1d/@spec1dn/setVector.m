function [ obj ] = setVector( obj,vec )
%SETVECTOR Summary of this function goes here
%   Detailed explanation goes here
    
validateattributes(vec,{'numeric'},{'size',[1 4]})
    
    if isfield(obj.resolution,'dqh')
        obj.resolution.dqh = vec;
    end
    if isfield(obj.resolution,'EXP')
        obj.resolution.ResCal.DQH = vec(1);
        obj.resolution.ResCal.DQK = vec(2);
        obj.resolution.ResCal.DQL = vec(3);
        obj.resolution.ResCal.DW  = vec(4);
    end
end

