function [varargout] = publishToGithub(file,varargin ) % post_name assets

p = inputParser;
p.addRequired('Filename');
p.addParameter('title',sprintf('Output-%i', round(1e5*rand(1)) ),@ischar);
p.addParameter('tags',[],@(x) iscell(x) || ischar(x));
p.addParameter('isdataset',false,@logical)
p.addParameter('execute',false,@logical)
p.addParameter('comments',true,@logical)
p.addParameter('description',[],@ischar)
p.addParameter('picture',[])
p.parse(file,varargin{:});

param = p.Results;

% keys = {'title','tags','description'};
%
% param = struct( 'title', sprintf('Output-%i', round(1e5*rand(1)) ) , ...
%     'tags', [],...
%     'isdataset', false, ...
%     'isreport', false, ...
%     'execute', false,...
%     'description',[]);
%
% for ii = 1 : numel( varargin )
%     if ischar( varargin{ii})
%         switch varargin{ii}
%             case 'title'
%                 param.title = varargin{ii+1};
%             case 'tags'
%                 param.tags = varargin{ii+1};
%             case 'description'
%                 param.description = varargin{ii+1};
%         end
%     end
% end
%
% %% Find the first time a paramter is called
%
% lastid = numel( varargin );
% if nargin > 1
%     cid = find( cellfun( @(x)ischar( x ), varargin ) );
%     if numel(cid) == 0
%         lastid = numel( varargin );
%     else
%         lastid = find( ismember( {varargin{cid}}, keys), 1,'first');
%         if numel(lastid) > 0
%             lastid = cid(lastid)-1;
%         else
%             lastid = numel( varargin );
%         end
%     end
% else
%     lastid = 1;
% end

persistent base_dir pid
if isempty(base_dir)
    base_dir = fullfile(sdext.getpref('libroot').val,'SpectraDocumentation',filesep);
end
if ~isempty(pid) % check for error in pid
    if isunix
        if ismac
            if system(sprintf('ps -fax | grep %i | grep -v grep',pid))
                pid = [];
            end
        else
            if system(sprintf('ps fax | grep %i | grep -v grep',pid))
                pid = [];
            end
        end
    else
        warning('Process information can not be obtained on a windows system.')
    end
else
    if isunix
        if ismac
            [code, temp]= system('ps -fax | grep jekyll | grep -v grep | awk ''{print $2}''');
        else
            [code, temp]= system('ps fax | grep jekyll | grep -v grep | awk ''{print $1}''');
        end
    else
        warning('Process information can not be obtained on a windows system.')
    end
    if ~isempty(temp)
        pid = str2double(temp);
    end
end
%% Post Type
% Determine the post type
% If varargin{1} is
%    a function or structure - then a dataset is generated
%    a script - then a report is generated
% This process selects the next steps in scripting the post.

if iscell( param.Filename ) || isstruct( param.Filename )
    % dataset
    param.isdataset = true;
elseif ischar( param.Filename )
    if strcmpi('start_server',param.Filename)
        code = system(sprintf('cd %s\n bundle exec jekyll serve &',base_dir));
        if code ~=0
            error()
        end
        if isunix
            if ismac
                [code, temp]= system('ps -fax | grep jekyll | grep -v grep | awk ''{print $2}''');
            else
                [code, temp]= system('ps fax | grep jekyll | grep -v grep | awk ''{print $1}''');
            end
        else
            warning('Process information can not be obtained on a windows system.')
        end
        if code ~=0
            error()
        end
        pid = str2double(pid);
        fprintf('Started jekyll server on http://127.0.0.1:4000 with pid %i\n',pid)
        return
    elseif strcmpi('stop_server',param.Filename)
        if isempty(pid)
            error('Server is not started')
        else
            code = system(sprintf('kill %i',pid));
            if code ~=0
                error()
            else
                pid = [];
            end
            fprintf('Killed jekyll server with pid %i\n',pid)
        end
        return
    elseif strcmpi('base_dir',param.Filename)
        base_dir = param.Filename;
        return
    else
        % report
        param.isreport = true;
        param.execute = true;
    end
    
else
    try
        f = functions( param.Filename );
        % dataset
        param.isdataset = true;
        param.execute = true;
    catch
        % error
    end
end

%% Parse Outputs

timenow = clock;
if param.isreport
    % Always default and send the contents to assets
    fmatpub = publish( param.Filename, 'outputDir', fullfile('.','assets') );
    param.layout = 'post';
    % elseif param.isdataset
    %     if param.execute
    %         [ varargout{1:nargout} ] =  param.Filename( varargin{2:lastid} );
    %
    %         varargout{1}.driver = func2str( varargin{1} );
    %     else
    %         varargout = { varargin{1} };
    %     end
    %     param.layout = 'dataset';
    % else
    %     error('something bad happened');
end

post_dir = fullfile(base_dir,'_posts');
pic_dir = fullfile(base_dir,'images');
pic_list = fullfile(pwd,'assets','*.png');
pic_names_l = dir(pic_list);
pic_names = {pic_names_l.name};
pic_names_old = pic_names;
if ~isdir(pic_dir)
    mkdir(pic_dir);
end
if ~isempty(pic_names)        
    copy_it = cellfun(@(x) copyfile(fullfile(pwd,'assets',x),fullfile(pic_dir, sprintf( '%04i-%02i-%02i-%s', timenow(1), timenow(2), timenow(3), regexprep(x,' ','-') ) )),{pic_names_l.name});
    newpid_loc = cellfun(@(x) fullfile(pic_dir, sprintf( '%04i-%02i-%02i-%s', timenow(1), timenow(2), timenow(3), regexprep(x,' ','-') ) ),{pic_names_l.name},'UniformOutput',0);
    pic_names = cellfun(@(x) sprintf( '%04i-%02i-%02i-%s', timenow(1), timenow(2), timenow(3), regexprep(x,' ','-') ),pic_names,'UniformOutput',0);
end
if ~isdir(post_dir)
    mkdir(post_dir);
end

%% Write pages
to_file = fullfile(post_dir, sprintf( '%04i-%02i-%02i-%s.html', timenow(1), timenow(2), timenow(3), regexprep( param.title,' ','-') ) );
fto = fopen( to_file ,'w');

%% YAML FRONT MATTER

fprintf( fto, '---\nlayout: %s\ntitle: "%s"\n', param.layout,regexprep( param.title,'-',' ') );

if numel( param.tags) > 0
    fprintf(fto,'tags:\n');
    fprintf(fto,'- %s\n', param.tags{:} );
end
if numel( param.description) > 0
    fprintf(fto,'description: "%s"\n', param.description );
end
if param.comments
    fprintf(fto,'comments: true\n');
end
if ~isempty(param.picture)
    fprintf(fto,'image:\n');
    fprintf(fto,'  feature: %s\n',param.picture);
end
%%

% close front matter header if is a report
if param.isreport
    fprintf(fto, '\n---\n');
    disp('Running Matlab''s publish function.');
    WebDat = fileread(fmatpub);
    
    cssig = MatlabCSS();
    
    for ii = 1 : numel(cssig)
        WebDat = regexprep( WebDat, cssig{ii}, '') ;
    end
    ind1 = strfind(WebDat,'html,body,div,span,applet');
    ind2 = strfind(WebDat,'table{border-collapse:collapse;border-spacing:0}')+length('table{border-collapse:collapse;border-spacing:0}');
    WebDat(ind1:ind2) = [];
    if ~isempty(pic_names)
        for ii = 1:length(pic_names)
            WebDat = regexprep(WebDat,sprintf('<img vspace="5" hspace="5" src="%s" alt="">',pic_names_old{ii}),sprintf('<figure><a href="{{ site.url }}/images/%1$s"><img src="{{ site.url }}/images/%1$s" alt="%2$s"></a></figure>',pic_names{ii},strrep(pic_names{ii},'.png','')));
            code = system(sprintf('cd %s\n git add %s',base_dir,newpid_loc{ii}));
            if code ~=0
                warning()
            end
        end
    end
    
    % Insert a javascript to reorganize image paths
    for ii = 1 : size(WebDat,1)
        fprintf( fto, '%s\n', WebDat(ii,:) );
    end
    % End report
    % else % STart dataset
    %     dskyfld = {'name','comment','image','url','link','description','include','html','driver'};
    %
    %     %% Variable names
    %     unique_variables = {};
    %     for ii = 1 : numel( varargout{1} )
    %         if isstruct( varargout{1} ) || isnumeric( varargout{1} )
    %             flds = fieldnames( varargout{1}(ii) );
    %             GetEl = @(x)varargout{1}(x);
    %         elseif iscell( varargout{1} );
    %             flds = fieldnames( varargout{1}{ii} );
    %             GetEl = @(x)varargout{1}{x};
    %         end
    %
    %         newfields = fieldnames(GetEl(ii));
    %         unreserved = logical(zeros( 1, numel( newfields ) ));
    %         for nn = 1 : numel( newfields )
    %
    %             if ~ismember( newfields{nn}, dskyfld )
    %                 unreserved(nn) = numel(getfield( GetEl(ii), newfields{nn} )) == 1;
    %             end
    %         end
    %
    %         flds = fieldnames(GetEl(ii));
    %         if any( unreserved )
    %             unique_variables = union( unique_variables, {flds{unreserved}});
    %         end
    %     end
    %
    %
    %     fprintf( fto, 'var:\n' );
    %     for ii = 1 : numel( unique_variables )
    %         fprintf( fto, '  - %s\n', unique_variables{ii} );
    %     end
    %
    %
    %     % Individual dataset level
    %     fprintf( fto,'data: \n' );
    %
    %     for ii = 1 : numel( varargout{1} )
    %
    %         if isstruct( varargout{1}(ii) ) || isnumeric( varargout{1}(ii) )
    %             flds = fieldnames( varargout{1}(ii) );
    %             GetEl = @(x)varargout{1}(x);
    %         elseif iscell( varargout{1}(ii) );
    %             flds = fieldnames( varargout{1}{ii} );
    %             GetEl = @(x)varargout{1}{x};
    %         end
    %
    %
    %
    %         for jj = 1 : numel( flds )
    %             fldstruct = getfield( GetEl(ii), flds{jj} );
    %             % All of this is performed on the first pass
    %             if jj == 1
    %                 % dataset metadata
    %                 initky = true;
    %                 for kys = 1 : numel(dskyfld)
    %
    %                     [ fldval cont ] = CheckGetField( GetEl(ii),  dskyfld{kys} );
    %                     if cont
    %
    %                         if initky
    %                             fprintf( fto, '- %s: \n', dskyfld{kys});
    %                             initky=false;
    %                         else
    %                             fprintf( fto, '  %s: \n', dskyfld{kys});
    %                         end
    %
    %                         if iscell( fldval )
    %                             for qq = 1 : numel( fldval )
    %                                 fprintf( fto, '  - %s \n', fldval{qq} );
    %                             end
    %                         else
    %                             fprintf( fto, '  - %s \n', fldval );
    %                         end
    %
    %                     end
    %                 end
    %                 if initky
    %                     fprintf( fto, '- metadata:\n');
    %                 else
    %                     fprintf( fto, '  metadata:\n');
    %                 end
    %             end
    %
    %
    %
    %             if ~ismember( flds{jj}, dskyfld )
    %                 fprintf( fto, '  - var: %s\n', flds{jj} );
    %                 if isnumeric( fldstruct ) && numel( fldstruct ) == 1
    %                     fprintf( fto, '    value: %f\n', fldstruct );
    %                 else
    %                     N = ndims( fldstruct );
    %                     fprintf( fto, '    dims: \n' );
    %                     for nn = 1 : N
    %                         fprintf( fto, '     - %i\n', size( fldstruct, nn) );
    %                     end
    %                     fprintf( fto, '    type: %s\n', class(fldstruct) );
    %                 end
    %
    %             end
    %         end
    %     end
    
end

%% CLOSE FILES

if param.isdataset
    fprintf(fto, '\n---\n');
end
fclose(fto);
code = system(sprintf('cd %s\n git add %s',base_dir,to_file));
if code ~=0
    warning()
end

if ~isempty(pid)
    pause(0.4)
    web('http://127.0.0.1:4000/')
end

end % END function

%% DEPENDENCIES

%% For Reports

function s = MatlabCSS()

s = {'html { min-height:100%; margin-bottom:1px; }';
    'html body { height:100%; margin:0px; font-family:Arial, Helvetica, sans-serif; font-size:10px; color:#000; line-height:140%; background:#fff none; overflow-y:scroll; }';
    'html body td { vertical-align:top; text-align:left; }';
    'h1 { padding:0px; margin:0px 0px 25px; font-family:Arial, Helvetica, sans-serif; font-size:1.5em; color:#d55000; line-height:100%; font-weight:normal; }';
    'h2 { padding:0px; margin:0px 0px 8px; font-family:Arial, Helvetica, sans-serif; font-size:1.2em; color:#000; font-weight:bold; line-height:140%; border-bottom:1px solid #d6d4d4; display:block; }';
    'h3 { padding:0px; margin:0px 0px 5px; font-family:Arial, Helvetica, sans-serif; font-size:1.1em; color:#000; font-weight:bold; line-height:140%; }';
    'a { color:#005fce; text-decoration:none; }';
    'a:hover { color:#005fce; text-decoration:underline; }';
    'a:visited { color:#004aa0; text-decoration:none; }';
    'p { padding:0px; margin:0px 0px 20px; }';
    'img { padding:0px; margin:0px 0px 20px; border:none; }';
    'p img, pre img, tt img, li img, h1 img, h2 img { margin-bottom:0px; }';
    'ul { padding:0px; margin:0px 0px 20px 23px; list-style:square; }';
    'ul li { padding:0px; margin:0px 0px 7px 0px; }';
    'ul li ul { padding:5px 0px 0px; margin:0px 0px 7px 23px; }';
    'ul li ol li { list-style:decimal; }';
    'ol { padding:0px; margin:0px 0px 20px 0px; list-style:decimal; }';
    'ol li { padding:0px; margin:0px 0px 7px 23px; list-style-type:decimal; }';
    'ol li ol { padding:5px 0px 0px; margin:0px 0px 7px 0px; }';
    'ol li ol li { list-style-type:lower-alpha; }';
    'ol li ul { padding-top:7px; }';
    'ol li ul li { list-style:square; }'};

s = cellfun( @(x)strtrim(x), cellstr(s), 'UniformOutput',false );

end


%% For Datasets
function [ value cont ] = CheckGetField( S, fld )
if isfield( S, fld)
    value = getfield( S, fld );
    cont = true;
else
    value = nan;
    cont = false;
end
end