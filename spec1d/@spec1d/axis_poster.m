function varargout=axis_poster(varargin)
opt={};
for i=1:length(varargin)
    temp=varargin{i};
    if isnumeric(temp)
        axes(temp)
    else
        opt{i}=temp;
    end
end
set(gca,'Fontname','Helvetica','Fontsize',21,'LineWidth',2);

if ~isempty(opt)
    set(gca,opt{:})
end
