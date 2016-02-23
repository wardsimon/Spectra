function varargout = plot(varargin)
%
% function plot(varargin)
%
% SPEC1D/PLOT Plot for spectra
%
% Usage: 1. Simple plot or overlay of spectra: >> plot(s1,s2,s3)
%                                              >> plot(s1,10.*s2)
%                                              >> plot(s1,10.*s2,'r-')
%                                              >> plot(s1,10.*s2,'r-','Marker','s')
%        2. Possibility to choose axes type:   >> plot(s1,'semilogy',1)
%                                              >> plot(s1,'semilogx',1)
%                                              >> plot(s1,'loglog',1)
%        !!! EXPERIMENTAL !!!
%        4. Create a log to go with the plot. This means we can trace
%           files and 'remember' how data was analysed.
%
% Simon Ward 27/01/2016

s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin(s_ind) = [];
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('s_in',@iscell);
p.addParamValue('semilogx',0,@isnumeric);
p.addParamValue('semilogy',0,@isnumeric);
p.addParamValue('loglog',0,@isnumeric);
p.addParamValue('semilog',0,@isnumeric);
p.addParamValue('trace',0,@isnumeric);
p.addParamValue('tLenght',NaN,@isnumeric);
% Try to deal with the 'ro' arguments etc...
% This is a total hack but matlab does not have an alternative.
if mod(length(varargin),2)
    p.addOptional('plot_spec',[],@ischar)
    p.parse(s,varargin{1},varargin{2:end});
    warning('off','MATLAB:structOnObject');
    p = struct(p);
    warning('on','MATLAB:structOnObject');
else
    p.parse(s,varargin{:});
    warning('off','MATLAB:structOnObject');
    p = struct(p);
    warning('on','MATLAB:structOnObject');
    p.Results.plot_spec = [];
end


s = [p.Results.s_in{:}];

plot_opt = p.Unmatched;
plot_spec = p.Results.plot_spec;

if p.Results.trace && sdext.getpref('experimental').val
    INFO = 'STACK TRACE\n-----------\n';
    
    [st, ~] = dbstack('-completenames');
    f = 0;
    for i = 1:length(st)
        [pathstr, ~, ~] = fileparts(st(i).file);
        % This assumes we are working from the home directory
        if ~isempty(strfind(pathstr,cast(java.lang.System.getProperty('user.home'),'char'))) && f == 0;
            rdir = pathstr;
            f = 1;
        end
        file_list{i} = st(i).file;
        line_list(i) = st(i).line;
    end
    [file_list, ind] = unique(file_list);
    INFO = sprintf('%sUser base directory is assumed to be: %s\nLine\tFile',INFO,rdir);
    for i = ind
        INFO = sprintf('%s\n%i\t%s\n',INFO,line_list(i),file_list{i});
    end
end

%% Plot the data
%  Work out limits on graph
min_x = min(arrayfun(@(x) min(x.x),struct(s)));
max_x = max(arrayfun(@(x) max(x.x),struct(s)));
min_y = min(arrayfun(@(x) min(x.y - x.e),struct(s)));
max_y = max(arrayfun(@(x) max(x.y + x.e),struct(s)));

if ishold(gcf) && ~isempty(get(gcf,'Children'))
    temp_x = get(gca,'XLim');
    temp_y = get(gca,'YLim');
    if temp_x(1) < min_x
        min_x = temp_x(1);
        pm(1) = 0.05*(max_x - min_x);
    else
        pm(1) = 0;
    end
    if temp_x(2) > max_x
        max_x = temp_x(2);
        pm(1) = 0.05*(max_x - min_x);
    end
    if temp_y(1) < min_y
        min_y = temp_y(1);
        pm(2) = 0.05*(max_y - min_y);
    else
        pm(2) = 0;
    end
    
    if temp_y(2) > max_y
        max_y = temp_y(2);
        pm(2) = 0.05*(max_y - min_y);
    end
    axis([min_x max_x min_y max_y])
    c_pos = length(findobj(get(gca,'Children'),'-regexp','Tag','[^'']'));
else
    if max(arrayfun(@(x) length(x.x),struct(s))) == 1
        pm = [0.05 0.05];
    else
        pm = [0.05*(max_x - min_x) 0.05*(max_y - min_y)];
    end
    c_pos = 0;
end

axis([min_x-pm(1) max_x+pm(1) min_y-pm(2) max_y+pm(2)])

hout  = [];
hbout = [];
hfout = [];

% Finally plot something
g = gcf;
held = ishold(g);
markerorder = {'o','square','diamond','^','>','<','p'};
colourorder = get(gca,'ColorOrder');
hold on

for i = 1:length(s)
    %----- Plot graph and error bars
    x = s(i).x(:); y = s(i).y(:); err = s(i).e(:); yfit = s(i).yfit(:);
    
    %----- build up nan-separated vector for bars
    ytop = (y + err)';
    ybot = (y - err)';
    if isnan(p.Results.tLenght)
        if length(x) == 1
            tee = 0.075;
        else
            tee = 0.075*min(diff(x));  % make tee distance for error bars
        end
    else
        if p.Results.tLenght > 0
            tee = p.Results.tLenght*min(diff(x)); % Allow for user scaling
        else
            tee = -1*p.Results.tLenght; % Allow for user supplied value
        end
    end
    xleft = (x-tee)';
    xright = (x+tee)';
    x_eb = cell2mat(arrayfun(@(xl, xr, x) [xl xr NaN x x NaN xl xr NaN],xleft,xright,x(:)','UniformOutput',0));
    y_eb = cell2mat(arrayfun(@(yb, yt) [yb yb NaN yb yt NaN yt yt NaN],ybot,ytop,'UniformOutput',0));
    %     for j = 1:length(x)
    %         x_eb((((j-1)*9)+1):(j*9)) = [xleft(j) xright(j) NaN x(j) x(j) NaN xleft(j) xright(j) NaN];
    %         y_eb((((j-1)*9)+1):(j*9)) = [ybot(j) ybot(j) NaN ybot(j) ytop(j) NaN ytop(j) ytop(j) NaN];
    %     end
    hle = plot(x_eb,y_eb,'Color',[128 128 128]/255,'LineStyle','-','LineWidth',1,'Marker','None');
    hEGroup = hggroup;
    set(hle,'Parent',hEGroup)
    
    %----- Do the actual plotting
    hll=plot(x,y,'MarkerSize',7,'MarkerEdgeColor','None',...
        'Tag',num2str(i),'DisplayName',sprintf('Dataset %i',i));
    de_opt = parse_opts(plot_spec);
    set(hll,de_opt{:})
    
    if ~isempty(plot_opt)
        f = fieldnames(plot_opt);
        for j = 1:length(f)
            set(hll,f{j},plot_opt.(f{j}))
        end
    end
    if ~isempty(yfit);
        hlf=plot(x,yfit,'Color',hsv2rgb([1 1 0.5].*rgb2hsv(get(hll,'MarkerFaceColor'))),...
            'LineStyle','-','LineWidth',2,'Marker','none','DisplayName',sprintf('Datafit %i',i));
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
    
    hout  = [hout;hll];
    hbout = [hbout;hle];
    hfout = [hfout;hlf];
end

if ~held
    hold off
end

% Create outputs
if nargout == 1
    % We only want points
    varargout{1} = hout;
elseif nargout == 2
    % We want points and error bars
    varargout{1} = hout;
    varargout{2} = hbout;
    % We want points, error bars and fit
elseif nargout == 3
    varargout{1} = hout;
    varargout{2} = hbout;
    varargout{3} = hfout;
end

set(gca,'Box','on')

%----- Make title

set(get(gca,'Title'),...
    'String',s(1).datafile,...
    'Tag','plot_text_title',...
    'FontName','Helvetica',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(get(gca,'Xlab'),...
    'String',s(1).x_label,...
    'Tag','plot_text_xlabel',...
    'FontName','Helvetica',...
    'Fontsize',14,...
    'buttondownfcn','edtext',...
    'handlevisibility','on',...
    'interpreter','none');

set(get(gca,'Ylab'),...
    'String',s(1).y_label,...
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
if any([p.Results.semilogx p.Results.semilogy p.Results.loglog])
    if p.Results.loglog
        p.Results.semilogx = 1;
        p.Results.semilogy = 1;
    end
    if p.Results.semilogx
        set(gca,'Xscale','log');
    end
    if p.Results.semilogy
        set(gca,'Yscale','log');
    end
end

if p.Results.trace && experimental
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

% Put axis on the top
set(gca,'Layer','Top')


    function default = parse_opts(in)
        ls = {'-','--',':','-.'};
        ma = {'o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'};
        co = {'y' 'm' 'c' 'r' 'g' 'b' 'w','k'};
        
        default = {'LineStyle','None','Marker',markerorder{1+mod(i-1,7)},'MarkerFaceColor',colourorder(1+mod(i-1 + c_pos,7),:)};
        
        for q = 1:length(in)
            c = in(q);
            if any(strcmp(c,ls))
                if length(in)>q
                    if any(strcmp(in(q+1),{'-','.'}))
                        c = in(q:(q+1));
                        in(q+1) = ' ';
                    end
                end
                default{2} = c;
                default{7} = 'Color';
                default{8} = default{6};
            end
            if any(strcmp(c,ma))
                default{4} = c;
            end
            if any(strcmp(c,co))
                default{6} = c;
                default{7} = 'Color';
                default{8} = c;
            end
        end
    end
end
