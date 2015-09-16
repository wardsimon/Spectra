function wait(varargin)
    % Improved pause function. 
    for i=1:length(varargin)
        for j=1:length(varargin{i})
            java.lang.Thread.sleep(varargin{i}(j)*1000);
        end
    end