function varargout = plot(varargin)
k=1;
name_given=0;
for i=1:nargin
    if isa(varargin{i},'ALPS')
        for j=1:length(varargin{i})
            ALPS(k)=varargin{i}(j);
            k=k+1;
        end
    else
        if isa(varargin{i},'cell')
            x_name=varargin{i}{1};
            y_name=varargin{i}{2};
            name_given=1;
        end
    end
end

plot_dwin

if (k-1)<=3
    colororder=[0 0 1; 1 0 0; 0 1 0];
else
    colororder=jet(k-1);
end

hold on
for i=1:length(ALPS)
    if name_given
        if strcmp(x_name,ALPS(i).xname)
            xn=ALPS(i).xname;
            if strcmp(y_name,ALPS(i).yname)
                p(i)=plot(ALPS(i).x,ALPS(i).y,'Color',colororder(i,:));
                yn=ALPS(i).yname;
            elseif strcmp(y_name,ALPS(i).zname)
                yn=ALPS(i).zname;
                p(i)=plot(ALPS(i).x,ALPS(i).z,'Color',colororder(i,:));
            end
        elseif strcmp(x_name,ALPS(i).yname)
            xn=ALPS(i).yname;
            if strcmp(y_name,ALPS(i).xname)
                p(i)=plot(ALPS(i).y,ALPS(i).x,'Color',colororder(i,:));
                yn=ALPS(i).xname;
            elseif strcmp(y_name,ALPS(i).zname)
                yn=ALPS(i).zname;
                p(i)=plot(ALPS(i).y,ALPS(i).z,'Color',colororder(i,:));
            end
        else
            xn=ALPS(i).zname;
            if strcmp(y_name,ALPS(i).xname)
                p(i)=plot(ALPS(i).z,ALPS(i).x,'Color',colororder(i,:));
                yn=ALPS(i).xname;
            elseif strcmp(y_name,ALPS(i).yname)
                p(i)=plot(ALPS(i).z,ALPS(i).y,'Color',colororder(i,:));
                yn=ALPS(i).yname;
            end
        end
    else
        p(i)=plot(ALPS(i).x,ALPS(i).y,'Color',colororder(i,:));
        xn=ALPS(i).xname;
        yn=ALPS(i).yname;
    end
end
hold off

%----- Make title

set(get(gca,'Title'),...
    'String',ALPS(end).observable,...
    'Tag','plot_text_title',...
    'FontName','Times',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(get(gca,'Xlab'),...
    'String',xn,...
    'Tag','plot_text_xlabel',...
    'FontName','Times',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(get(gca,'Ylab'),...
    'String',yn,...
    'Tag','plot_text_ylabel',...
    'FontName','Times',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(gca,'FontName','Times','FontSize',14);
set(gca,'Xgrid','on','Ygrid','on');
set(gca,'ButtonDownFcn','scribeaxesdlg(gca);');

if nargout==1
    varargout=p;
else
    varargout=[];
end

%======================================================================================
function plot_dwin
%
% MATLAB function to draw data window and add menus
%
%
%

% HMR 2002/1
%hplot=findobj('Tag','s1d_DataWindow');
hplot=gcf;

if ~isempty(hplot)
    figure(hplot)
    if strcmp(get(gca,'nextplot'),'replace')
        cla
    end
else
    hplot=figure('Position',[50 50 500 400],...
        'Tag','s1d_DataWindow',...
        'Visible','on',...
        'Name','Spectra Plot',...
        'Windowbuttondownfcn','edtext(''hide'')',...
        'Menubar','None',...
        'NumberTitle','off');
    
    %---------- Print menu-----------------------------
    
    hprt=uimenu(hplot,'Label','Print');
    uimenu(hprt,'Label','Print figure',...
        'Callback',['hplot=findobj(''Tag'',''s1d_DataWindow'');'...
        'figure(hplot);'...
        'set(gca,''units'',''points'');'...
        'pos=get(gca,''Position'');'...
        'set(gcf,''papertype'',''a4letter'');'...
        'set(gcf,''units'',''centimeters'');'...
        'set(gcf,''paperunits'',''centimeters'');'...
        'p=get(gcf,''position'');'...
        'width=p(3); height=p(4);'...
        'paper=get(gcf,''papersize'');'...
        'ppos=[(paper(1)-width)/2 (paper(2)-height)/1.25 width height];'...
        'set(gcf,''paperposition'',ppos);'...
        'printdlg(gcf);']);
    
    uimenu(hprt,'Label','Save figure',...
        'Callback','[f,p] = uiputfile(''*.m'',''Save Figure As...'');if f,eval([''print -dmfile '' p f ]); disp(''Figure saved''); end');
    
end
