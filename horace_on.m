function path=horace_on(non_default_path)
%  safely switches horace on
%  horace_on()                         -- calls horace with default settings
%  horace_on(non_default_horace_path)  -- calls horace with non-default horace folder;
%
%
% $Revision: 1733 $ ($Date: 2010-06-23 17:57:14 +0100 (Wed, 23 Jun 2010) $) 
%
%	default_horace_path = [matlabroot() filesep 'toolbox' filesep 'Horace_SVN'];


try
    libroot = sdext.getpref('libroot').val;
catch
    libroot = fullfile(matlabroot(),'toolbox','mtools');
end

default_herbert_path = fullfile(libroot,'Herbert');
default_horace_path = fullfile(libroot,'Horace');

prefer_herbert=true;
%
if exist('non_default_path','var') && (strcmpi(non_default_path,'where') || strcmpi(non_default_path,'which'))
    path = find_default_path(default_horace_path);   
    return;
end

warn_state=warning('off','all');    % turn of warnings (so dont get errors if remove non-existent paths)
try
   horace_off();
catch
end
warning(warn_state);    % return warnings to initial state
if ~isempty(default_herbert_path)
	herbert_on(default_herbert_path);
end

libisis_initated=~isempty(which('libisis_init.m'));
herbert_initated=~isempty(which('herbert_init.m'));

% if neither libisis not herbert are initated, try to init something. 
if ~(libisis_initated||herbert_initated)
    if prefer_herbert
        try
            herbert_on();               
        catch
            try
                libisis_on();
            catch
            end
            
        end
    else
        try
            libisis_on();
        catch
            try
                herbert_on();               
            catch
            end
        end
    end
end

% init horace 
if nargin==1 
	start_app(non_default_path);	
else
	start_app(default_horace_path);
end
path = fileparts(which('horace_init.m'));


function start_app(path)
addpath(path);
horace_init;

function path =find_default_path(her_default_path)
path = which('horace_init.m');
if isempty(path)
    path = her_default_path;
    if ~exist(fullfile(path,'horace_init.m'),'file')
        path='';
    end
else
    path=fileparts(path);
end


 
function path=herbert_on(non_default_path)
% The function intended to swich herbert on and
% return the path were herbert is resided or 
% empty string if herbert has not been found
%
%
% The function has to be present in Matlab search path 
% and modified for each machine to know default herbert location
%
%Usage:
%>>path=herbert_on(); 
%       enables herbert and initiates herbert default search path
%>>path=herbert_on('where'); 
%       reports current location of herbert or empty if not found
%
%>>path=herbert_on('a path'); 
%       initiates herbert on non-default search path
%
%
%
%
her_default_path=[fileparts(which('horace_on.m')) filesep 'Herbert'];
%
if exist('non_default_path','var') && (strcmpi(non_default_path,'where') || strcmpi(non_default_path,'which'))
    path = find_her_default_path(her_default_path);   
    return;
end
if nargin==1 
    start_herbert(non_default_path);    
else
    start_herbert(her_default_path);    
end
path = fileparts(which('herbert_init.m'));


function start_herbert(path)

try
    herbert_off;
catch
end
addpath(path);
herbert_init;

function path =find_her_default_path(her_default_path)
path = which('herbert_init.m');
if isempty(path)
    path = her_default_path;
    if ~exist(fullfile(path,'herbert_init.m'),'file')
        path='';
    end
else
    path=fileparts(path);
end
