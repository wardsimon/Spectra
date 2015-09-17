%% SYSTEM ILL STARTUP FILE
% File        : startup.m
% Description : Matlab libraries startup file
% Author      : Simon Ward
% Date        : 10/04/2013
% Changes     : Added Unix/Windows compatibilty
%               try/catch renaming
%               Use java as it is blind to your system.
%               Robust path additions


% Use Java, it's more reliable!
st_home = cast(java.lang.System.getProperty('user.home'),'char');
st_user = cast(java.lang.System.getProperty('user.name'),'char');
try
    st_comp = cast(java.net.InetAddress.getLocalHost.getHostName,'char');
catch
    st_comp = 'unavailable';
end
d = [st_home filesep 'Documents' filesep 'MATLAB'];
if isdir(d)
    cd(d)
else
    d = [getenv('HOME') filesep 'MATLAB'];
    if isdir(d)
        cd(d)
    end
end

% Start a diary
if exist([d filesep 'matlab.log'],'file') == 2
    try
        if exist(sprintf('matlab_%s.log',date),'file') == 2
            if ispc
                system(sprintf('matlab.log >> matlab_%s.log',date));
            elseif isunix
                system(sprintf('cat matlab.log >> matlab_%s.log',date));
            end
        else
            copyfile('matlab.log', sprintf('matlab_%s.log',date))
        end
        delete('matlab.log');
        fprintf('Previous Matlab Log file is renamed as matlab_%s.log\n',date)
    catch ME
        fprintf('Renaming of previous Matlab log file has failed: \n%s\n',ME.message)
    end
end
diary([d filesep 'matlab.log'])
disp([ 'matlab.log started on : ' datestr(now) ])

% THIS FILE CONTAINS the libroot!
if exist('startuser.m','file') == 2
    eval('startuser');
end

% Add path to libraries
% ! NOTE ! We do not have to add files to path if this is correct!
% Go to default place ???
try
    libroot = getpref('mtools','libroot');
    if ~exist('libroot',1)
        libroot = uigetdir(st_home,'Select mtools root directory');
        if libroot == 0
            error('Without an mtools directory, here be dragons!')
        else
            setpref('mtools','libroot',libroot)
            choice = questdlg('Would you to enable experimental features?', ...
                'Enable Extras', ...
                'Yes','No','No');
            % Handle response
            switch choice
                case 'No'
                    setpref('mtools','experimental',0)
                case 'Yes'
                    setpref('mtools','experimental',1)
            end
        end
    end
catch
    libroot = uigetdir(st_home,'Select mtools root directory');
    if libroot == 0
        error('Without an mtools directory, here be dragons!')
    else
        setpref('mtools','libroot',libroot)
        choice = questdlg('Would you to enable experimental features?', ...
            'Enable Extras', ...
            'Yes','No','No');
        % Handle response
        switch choice
            case 'No'
                setpref('mtools','experimental',0)
            case 'Yes'
                setpref('mtools','experimental',1)
        end
    end
end

addpath(genpath(fullfile(libroot,'Spectra')))

if ~isempty(st_user) && ~isempty(st_comp)
    fprintf( 'Welcome %s@%s! ILL files are loaded\n',st_user,st_comp);
end

% If you have the curve fitting toolbox horace and herbert overwrite some functions
v=ver;
if any(strcmp('Curve Fitting Toolbox', {v.Name}))
    try
        herbert_off
        horace_off
        fprintf('!! Horace and Herbert have been disabled !!\n!!  Use horace_on, herbert_on to enable  !!\n')
    end
end

% Clear up anything I have created
clear all
close all

% Load previous workspace
if exist('matlab.mat','file') == 2
    disp('Getting previous workspace data from matlab.mat')
    load matlab.mat
end

% Set defaults
set(0,'DefaultFigurePaperUnits','centimeters');
set(0,'DefaultFigurePaperType','A4');
