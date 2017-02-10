function OK = checkGpuMemory(s)
% Calculates the free memory on the GPU for a spec1d object
%
% OK = CHECKGPUMEMORY(obj)
%
% Input:
%
% obj   spec1d class object
%
% Output:
% OK    boolean value
%
% This function checks if there is enough space on the GPU for obj. Returns
% true if it can be pushed, false otherwise.
%

try
    d = gpuDevice; % Try to get GPU device.
    memA = d.AvailableMemory;  % Check available memory
    memD = 3*length(s.x)*8; % We will use this much memory
    
    % Check if enough free space...
    if memD < memA
        OK = 1;
    else
        OK = 0;
    end
catch
    OK = 0;
end

end