function OK = checkGpuMemory(s)
try
    d = gpuDevice;
    memA = d.AvailableMemory;  % Check available memory
    memD = 3*length(s.x)*8; % We will use this much memory
    if memD < memA
        OK = 1;
    else
        OK = 0;
    end
catch
    OK = 0;
end
end