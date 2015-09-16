function varargout=ylabel_poster(varargin)

st=0;
opt={};
for i=1:length(varargin)
    temp=varargin{i};
    if isnumeric(temp)
        axes(temp)
    elseif ischar(temp) & ~st
        xla=temp;
        st=1;
    else
        opt{i}=temp;
    end
end

varargout{1}=ylabel(xla,'Fontname','Helvetica','Fontsize',21,'Interpreter','tex');
if ~isempty(opt)
    set(varargout{1},opt{:})
end
if nargout==0
    varargout{1}=[];
end