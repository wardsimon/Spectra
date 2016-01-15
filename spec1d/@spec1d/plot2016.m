function varargout = plot2016(varargin)
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
% Simon Ward 15/01/2016

p = inputParser;
p.addRequired('s_in',@(x) isa(x,'spec1d'));
p.addOptional('add_s',[],@(x) isa(x,'spec1d'))
p.addParameter('semilogx',0,@isnumeric);
p.addParameter('semilogy',0,@isnumeric);
p.addParameter('loglog',0,@isnumeric);
p.addParameter('semilog',0,@isnumeric);
p.addParameter('trace',0,@isnumeric);
p.addParameter('tLenght',NaN,@isnumeric);

p.parse(varargin{:});

s = [p.Results.s_in(:) p.Results.add_s(:)];
plot_opt = p.Unmatched;

if p.Results.trace && sdex.getpref('experimental').val
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

%% Plot the data
%  Work out limits on graph
min_x = Inf;
max_x = -Inf;
min_y = Inf;
max_y = -Inf;

for i = 1:length(s)
    if min(s(i).x) < min_x; min_x = min(s(i).x); end
    if max(s(i).x) > max_x; max_x = max(s(i).x); end
    if min(s(i).y - s(i).e) < min_y; min_y = min(s(i).y - s(i).e); end
    if max(s(i).y + s(i).e) > max_y; max_y = max(s(i).y + s(i).e); end
end
pm = [0.05*(max_x - min_x) 0.05*(max_y - min_y)];
axis([min_x-pm(1) max_x+pm(1) min_y-pm(2) max_y+pm(2)])


hout=[];
hbout=[];
hfout=[];

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
        tee = 0.075*min(diff(x));  % make tee distance for error bars
    else
        if p.Results.tLenght > 0
            tee = p.Results.tLenght*min(diff(x)); % Allow for user scaling 
        else
            tee = -1*p.Results.tLenght; % Allow for user supplied value
        end
    end
    xleft = (x-tee)';
    xright = (x+tee)';
    x_eb = zeros(9*length(x),1);
    y_eb = x_eb;
    for j = 1:length(x)
        x_eb((((j-1)*9)+1):(j*9)) = [xleft(j) xright(j) NaN x(j) x(j) NaN xleft(j) xright(j) NaN];
        y_eb((((j-1)*9)+1):(j*9)) = [ybot(j) ybot(j) NaN ybot(j) ytop(j) NaN ytop(j) ytop(j) NaN];
    end
    hle = plot(x_eb,y_eb,'Color',[128 128 128]/255,'LineStyle','-','LineWidth',1,'Marker','None');
    hEGroup = hggroup;
    set(hle,'Parent',hEGroup)
    
    %----- Do the actual plotting
    hll=plot(x,y,'LineStyle','None','Marker',markerorder{1+mod(i-1,7)},'MarkerFaceColor',colourorder(1+mod(i-1,7),:),...
        'MarkerSize',7,'MarkerEdgeColor','None','Tag',num2str(i),'DisplayName',sprintf('Dataset %i',i));
    
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

