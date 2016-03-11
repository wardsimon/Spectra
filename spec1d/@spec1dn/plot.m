function varargout = plot(s,varargin)

% Plot using superclass and get all plot results
varargout = cell(3,1);
[varargout{:}]  = plot@spec1d(s,varargin{:});

% Attach ID's to the plot points.
for i = 1:length(varargout{1})
   varargout{1}(i).UserData.ID = s(i).ident; 
end
% Attach listner to the points
set(varargout{1},'ButtonDownFcn',@(ob, event) resolutionMouseButtonPress(s,ob,event));

% Set return values.
if nargout >=1
    varargout = varargout{1:nargout};
else
    varargout = {};
end

end