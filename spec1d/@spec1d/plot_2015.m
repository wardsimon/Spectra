function varargout=plot(varargin)
%
% function plot(varargin)
%
% SPEC1D/PLOT Plot for spectra
%
% Usage: 1. Simple plot or overlay of spectra: >> plot(s1,s2,s3)
%                                              >> plot(s1,10.*s2)
%        2. Possibility to choose axes type:   >> plot(s1,'semilogy')
%                                              >> plot(s1,'semilogx')
%                                              >> plot(s1,'loglog')
%        3. More plot options                  >> plot(s1,'r-*') etc.
%
%        !!! EXPERIMENTAL !!!
%        4. Create a log to go with the plot. This means we can trace
%           files and 'remember' how data was analysed.
%
% DFM 2.10.98
% HMR 21.5.99 Outputs handle to markers and to errorbars (if asked for).
% SW  <<.4.13 Re-written, copies matlab plot and now has nice figures


%% !!! SIMONS EXPERIMENTAL ! ! !
experimental = sdext.getpref('experimental').val;

% Logging features. Still under development but works
global logplot
if ~isempty(logplot)
    if logplot==1 && experimental
        INFO='STACK TRACE\n-----------\n';
        
        [st, ~]=dbstack('-completenames');
        f=0;
        for i=1:length(st)
            [pathstr, ~, ~] = fileparts(st(i).file);
            % This assumes we are working from the home directory
            if ~isempty(strfind(pathstr,cast(java.lang.System.getProperty('user.home'),'char'))) && f==0;
                rdir=pathstr;
                f=1;
            end
            file_list{i}=st(i).file;
            line_list(i)=st(i).line;
        end
        [file_list, ind]=unique(file_list);
        INFO=sprintf('%sUser base directory is assumed to be: %s\nLine\tFile',INFO,rdir);
        for i=ind
            INFO=sprintf('%s\n%i\t%s\n',INFO,line_list(i),file_list{i});
        end
    end
end

%% Sort out arguments
%----- Draw data window
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot
%          c     cyan          +     plus               --    dashed
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%          w     white         v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram
plot_dwin;

oldvarargin=varargin;
varargin={}; s={};
plot_options={};
use_log=0;
color=[]; marker=[]; lines=[]; ml_pt=1; not_skip=1;
for i=1:length(oldvarargin)
    if isa(oldvarargin{i},'spec1d')
        for n=1:length(oldvarargin{i})
            s{end+1}=oldvarargin{i}(n);
        end
    elseif ischar(oldvarargin{i})
        % Check for log axis
        if sum(strcmp(oldvarargin{i},{'semilogy','semilogx','loglog'}))>0 && use_log==0
            log_string=oldvarargin{i};
            use_log=1;
            % check for Matlab plot options
        elseif length(oldvarargin{i})<=3 && ml_pt
            for j=1:length(oldvarargin{i})
                switch lower(oldvarargin{i}(j))
                    case 'b', color  = [0 0 1];
                    case 'g', color  = [0 1 0];
                    case 'r', color  = [1 0 0];
                    case 'c', color  = [0 1 1];
                    case 'm', color  = [1 0 1];
                    case 'y', color  = [1 1 0];
                    case 'k', color  = [0 0 0];
                    case 'w', color  = [1 1 1];
                    case '.', marker = oldvarargin{i}(j);
                    case 'o', marker = oldvarargin{i}(j);
                    case 'x', marker = oldvarargin{i}(j);
                    case '+', marker = oldvarargin{i}(j);
                    case '*', marker = oldvarargin{i}(j);
                    case 's', marker = oldvarargin{i}(j);
                    case 'd', marker = oldvarargin{i}(j);
                    case 'v', marker = oldvarargin{i}(j);
                    case '^', marker = oldvarargin{i}(j);
                    case '<', marker = oldvarargin{i}(j);
                    case '>', marker = oldvarargin{i}(j);
                    case 'p', marker = oldvarargin{i}(j);
                    case 'h', marker = oldvarargin{i}(j);
                    case '-', lines  = oldvarargin{i}(j);
                    case ':', lines  = oldvarargin{i}(j);
                end
            end
            ml_pt=0;
            % Check for direct passing options
        elseif not_skip
            plot_options{end+1,1}=oldvarargin{i};
            try
                plot_options{end,2}=oldvarargin{i+1};
            catch
                error('Argument needed for option: %s',oldvarargin{i})
            end
            not_skip=0;
        end
        varargin{end+1}=oldvarargin{i};
    end
end


%% Plot the data
%  Work out limits on graph
min_x = Inf;
max_x = -Inf;
min_y = Inf;
max_y = -Inf;

for i = 1:length(s)
    if min(s{i}.x) < min_x; min_x = min(s{i}.x); end
    if max(s{i}.x) > max_x; max_x = max(s{i}.x); end
    if min(s{i}.y - s{i}.e) < min_y; min_y = min(s{i}.y - s{i}.e); end
    if max(s{i}.y + s{i}.e) > max_y; max_y = max(s{i}.y + s{i}.e); end
end
pm = [0.05*(max_x - min_x) 0.05*(max_y - min_y)];
axis([min_x-pm(1) max_x+pm(1) min_y-pm(2) max_y+pm(2)])


% Make the arguments the same form as the data
if isempty(color)
    if experimental
        colororder=colorwheel('steps',length(s),'s',0.75,'v',0.975);
    else
        if length(s)<=3
            colororder=[0 0 1; 1 0 0; 0 1 0];
        else
            colororder=jet(length(s));
        end
    end
else
    colororder=hsv2rgb([ones(length(s),1) linspace(0.3,1,length(s))' ones(length(s),1)].*repmat(rgb2hsv(color),length(s),1));
end

if isempty(marker)
    markerorder=str2mat('o','square','diamond','^','>','<','p');
else
    markerorder=repmat(marker,7,1);
end

if isempty(lines)
    join_points=0;
    lineorder=strvcat('-','--',':','-.');
else
    join_points=1;
    lineorder=repmat(lines,4,1);
end

hout=[];
hbout=[];
hfout=[];

% Finally plot something
for i=1:length(s)
    %----- Plot graph and error bars
    x=s{i}.x; y=s{i}.y; err=s{i}.e; yfit=s{i}.yfit;
    x=x(:); y=y(:); err=err(:); yfit=yfit(:);
    
    %----- build up nan-separated vector for bars
    ytop = (y + err)';
    ybot = (y - err)';
    tee = 0.075*min(diff(x));  % make tee distance for error bars
    xleft = (x-tee)';
    xright = (x+tee)';
    x_eb = zeros(9*length(x),1);
    y_eb = x_eb;
    for j = 1:length(x)
        x_eb((((j-1)*9)+1):(j*9)) = [xleft(j) xright(j) NaN x(j) x(j) NaN xleft(j) xright(j) NaN];
        y_eb((((j-1)*9)+1):(j*9)) = [ybot(j) ybot(j) NaN ybot(j) ytop(j) NaN ytop(j) ytop(j) NaN];
    end
    hle = plot(x_eb,y_eb,'Color',[128 128 128]/255,'LineStyle','-','LineWidth',1,'Marker','None');
    hold on
    hEGroup = hggroup;
    set(hle,'Parent',hEGroup)

    hll=plot(x,y,'color',colororder(i,:),...
        'LineStyle','None','Marker',deblank(markerorder(1+mod(i-1,7),:)),...
        'MarkerSize',7,'MarkerFaceColor',colororder(i,:),'MarkerEdgeColor','None','Tag',num2str(i),'DisplayName',sprintf('Dataset %i',i));
    
    if experimental
        set(hll,'ButtonDownFcn',@(ob, event) buttonTest(ob,event));
        text('Color',hsv2rgb([1 1 0.7].*rgb2hsv(colororder(i,:))), 'VerticalAlign', 'Middle','Visible','Off','Tag',num2str(i),'Interpreter','Tex');
        % Performance problems with high density data
        if length(x) < 250
%             set(gcf,'WindowButtonMotionFcn', @hoverCallback);
        end
    else
        set(hll,'ButtonDownFcn','editline(gco);');
    end
    
    
    if join_points
        set(hll,'LineStyle',deblank(lineorder(1+mod(i-1,4),:)))
    end
    if ~isempty(plot_options)
        for j=1:size(plot_options,1)
            def=plot_options{j,1};
            var=plot_options{j,2};
            % Check for color options
            if ~isempty(strfind(def,'Color')) || ~isempty(strfind(def,'color'))
                if ischar(var)
                    var=lower(var);
                    switch var
                        case {'b', 'blue'   }, var  = [0 0 1];
                        case {'g', 'green'  }, var  = [0 1 0];
                        case {'r', 'red'    }, var  = [1 0 0];
                        case {'c', 'cyan'   }, var  = [0 1 1];
                        case {'m', 'magenta'}, var  = [1 0 1];
                        case {'y', 'yellow' }, var  = [1 1 0];
                        case {'k', 'black'  }, var  = [0 0 0];
                        case {'w', 'white'  }, var  = [1 1 1];
                        otherwise
                            warning('Color options must be r,g,b,c,m,y,k,w. Not %s',var)
                            var=[0 0 1];
                    end
                end
            end
            t=linspace(1,0.1,length(s));
            var=hsv2rgb([1 t(i) 1].*rgb2hsv(var));
            set(hll,def,var);
        end
    end
    
    if ~isempty(yfit);
        hlf=plot(x,yfit,'color',...
            hsv2rgb([1 1 0.5].*rgb2hsv(colororder(i,:))),...
            'LineStyle','-','LineWidth',2,'Marker','none','DisplayName',sprintf('Datafit %i',i));
        if experimental
            text('Color',hsv2rgb([1 1 0.7].*rgb2hsv(colororder(i,:))),...
                'VerticalAlign', 'Middle','Visible','Off','Tag',num2str(i+length(s)),'Interpreter','Tex');
            set(hlf,'Tag',num2str(i+length(s)))
%             set(gcf,'WindowButtonMotionFcn', @hoverCallback);
        end
        set(hlf,'ButtonDownFcn','editline(gco);');
        set(hlf,'Parent',hEGroup)
    else
        hlf=[];
    end
    
    % Remove error and fit from legend
    if ishg2
        % HG2 Way
        hle.LegendDisplay = 'off';
    else
        try 
            % The undocumented way
            hasbehavior(hle,'legend',0)
        catch
        % Foolproof way!
        set(get(get(hEGroup,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        end
    end
    
    hout=[hout;hll];
    hbout=[hbout;hle];
    hfout=[hfout;hlf];
end
    

% Create outputs
if nargout==1
    % We only want points
    varargout{1}=hout;
elseif nargout==2
    % We want points and error bars
    varargout{1}=hout;
    varargout{2}=hbout;
    % We want points, error bars and fit
elseif nargout==3
    varargout{1}=hout;
    varargout{2}=hbout;
    varargout{3}=hfout;
end

set(gca,'Box','on')

%----- Make title

set(get(gca,'Title'),...
    'String',s{1}.datafile,...
    'Tag','plot_text_title',...
    'FontName','Helvetica',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(get(gca,'Xlab'),...
    'String',s{1}.x_label,...
    'Tag','plot_text_xlabel',...
    'FontName','Helvetica',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(get(gca,'Ylab'),...
    'String',s{1}.y_label,...
    'Tag','plot_text_ylabel',...
    'FontName','Helvetica',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','tex');

set(gca,'FontName','Helvetica','FontSize',14);

if ismac
    % For some reason this looks crap on a mac
    set(gca,'Xgrid','on','Ygrid','on','Box','On','LineWidth',1);
else
    set(gca,'Xgrid','on','Ygrid','on','Box','On','LineWidth',2);
end

%----- Change to log axes if asked
if use_log
    switch log_string
        case 'semilogy', set(gca,'Yscale','log');
        case 'semilogx', set(gca,'Xscale','log');
        case 'loglog'  , set(gca,'Xscale','log','Yscale','log');
        otherwise, disp('Not a valid axis type')
    end
end

if ~isempty(logplot)
    if logplot==1 && experimental
        file_list{end+1}=[tempname '.jpg'];
        print(gcf,'-djpeg',file_list{end})
        while exist(file_list{end},'file')~=2,pause(0.05),end
        file_list{end+1}=strrep(file_list{end},'jpg','mat');
        evalin('base', ['save(''' file_list{end} ''')']);
        while exist(file_list{end},'file')~=2,pause(0.05),end
        zip_name=[rdir filesep 'Archive_' strrep(sprintf('%0.1f',clock),'.','_') '.zip'];
        INFO=sprintf('%sSaved archive: %s\nWorkspace, files and figures are included.',INFO,zip_name);
        info_name=[tempdir filesep 'README.txt'];
        fid=fopen(info_name,'w');
        fprintf(fid,INFO);
        while exist(info_name,'file')~=2,pause(0.02),end
        fclose(fid);
        file_list{end+1}=info_name;
        zip(zip_name,file_list)
        while exist(zip_name,'file')~=2,pause(0.05),end
        s=get(gcf,'UserData');
        fid=fopen(zip_name);
        s.GenBinaryZip=fread(fid);
        fclose(fid);
        set(gcf,'UserData',s);
    end
end

% Put axis on the top
set(gca,'Layer','Top')
end
%----- Start of subroutines

%======================================================================================
function plot_dwin
%
% MATLAB function to draw data window and add menus
%
%
%
c_win=get(0,'CurrentFigure');
if isempty(c_win)
    c_win=figure;
end
s=struct('Figure',c_win);
set(c_win,'UserData',s);
end

function buttonTest(ob,event)
c_fig=findobj('Tag','SepData');
if isempty(c_fig)
    c_fig=figure('Tag','SepData');
    s1=subplot(2,2,[1 3]);
    x=get(ob,'XData');
    y=get(ob,'YData');
    p=plot(x,y);
    hold on
    set(p,'Color',get(ob,'Color'),'MarkerFaceColor',get(ob,'MarkerFaceColor'),'MarkerEdgeColor',get(ob,'MarkerEdgeColor'),'Marker',get(ob,'Marker'),'LineStyle',get(ob,'LineStyle'))
    axis([min(x)-(max(x)-min(x))*0.1 max(x)+(max(x)-min(x))*0.1 min(y)-(max(y)-min(y))*0.1 max(y)+(max(y)-min(y))*0.1]);
    s2=subplot(2,2,[2 4]);
    set(s2,'Visible','Off')
    t = uitable(...
        'Units','Normalized',...
        'BackgroundColor',[1 1 1;0.96078431372549 0.96078431372549 0.96078431372549],...
        'ColumnFormat',{  [] [] },...
        'ColumnEditable',[false false],...
        'ColumnName',{  'X'; 'Y' },...
        'ColumnWidth',{  'auto' 'auto' },...
        'Data',[x(:) y(:)],...
        'Position',get(s2,'Position'),...
        'UserData',[]);
    s=struct('FigureAxis',s1,'DataTable',t,'Parent',gcf);
    set(c_fig,'UserData',s)
else
    s=get(c_fig,'UserData');
    if isstruct(s)
        if all([isfield(s,'FigureAxis') isfield(s,'DataTable')])
            ax1=s.FigureAxis;
            x=get(ob,'XData');
            y=get(ob,'YData');
            p=plot(ax1,x,y);
            hold on
            set(p,'Color',get(ob,'Color'),'MarkerFaceColor',get(ob,'MarkerFaceColor'),'MarkerEdgeColor',get(ob,'MarkerEdgeColor'),'Marker',get(ob,'Marker'),'LineStyle',get(ob,'LineStyle'));
            old_data=get(s.DataTable,'Data');
            if size(old_data,1)>length(x)
                new_data=zeros(size(old_data,1),size(old_data,2)+2);
                new_data(:,1:size(old_data,2))=old_data;
                new_data(1:length(x),(end-1):end)=[x(:) y(:)];
            elseif size(old_data,1)<length(x)
                new_data=zeros(length(x),size(old_data,2)+2);
                new_data(1:size(old_data,1),1:size(old_data,2))=old_data;
                new_data(:,(end-1):end)=[x(:) y(:)];
            else
                new_data=[old_data x(:) y(:)];
            end
            
            set(s.DataTable,'Data',new_data,...
                'ColumnFormat',cell(1,length(old_data)+2),...
                'ColumnEditable',false(1,length(old_data)+2),...
                'ColumnName',{ repmat(['X';'Y'],size(old_data,2)/2 +1,1)})
        end
    end
end
end

function hoverCallback(src, evt)
% Grab the x & y axes coordinate where the mouse is
mousePoint = get(gca, 'CurrentPoint');
mouseX = mousePoint(1,1);
mouseY = mousePoint(1,2);
% Compare where data and current mouse point to find the data point
% which is closest to the mouse point
p_handles = findobj(gca,'Type','line','-not','Tag','');
t_handles = findobj(gca,'Type','text','-not','Tag','');
% If the distance is less than some threshold, set the text
% object's string to show the data at that point.
xl = get(gca, 'Xlim');
xrange = range(xl);
%     yrange = range(get(gca, 'Ylim'));
for i = 1:length(p_handles)
    for j = 1:length(t_handles)
        if str2double(get(t_handles(j),'Tag'))==str2double(get(p_handles(i),'Tag'))
            textHdl = t_handles(j);
            lineHdl = p_handles(i);
        end
    end
    try
        y_d = get(lineHdl,'Ydata');
        x_d = get(lineHdl,'Xdata');
        distancesToMouse = hypot(x_d - mouseX, y_d - mouseY);
        [val, ind] = min(abs(distancesToMouse));
        if abs(mouseX - x_d(ind)) < 0.05*xrange
            set(textHdl, 'String', sprintf('\\leftarrow {\\bf %0.5f} | %s',y_d(ind),get(lineHdl,'DisplayName')), 'Position', [x_d(ind), y_d(ind)])
            temp = get(textHdl,'Extent');
            if x_d(ind)+temp(3) > xl(2)
                set(textHdl,'HorizontalAlignment','Right','Visible','On')
                set(textHdl, 'String', sprintf('%s | {\\bf %0.5f} \\rightarrow ',get(lineHdl,'DisplayName'),y_d(ind)), 'Position', [x_d(ind), y_d(ind)])
            else
                set(textHdl,'HorizontalAlignment','Left','Visible','On')
            end
            %                 nd_handles = p_handles(p_handles ~= lineHdl);
            %                 for j = 1:length(nd_handles)
            %                    set(nd_handles(j),'Color',hsv2rgb([1 0.5 1].*rgb2hsv(get(nd_handles(j),'Color'))))
            %                 end
        else
            set(textHdl,'HorizontalAlignment','Left','Visible','Off')
        end
    catch
    end
end
end