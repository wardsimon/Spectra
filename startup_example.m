%% SPECTRA STARTUP FILE
% File        : startup.m
% Description : Matlab libraries startup file
% Author      : Simon Ward
% Date        : 13/01/2016
% Changes     : Added Unix/Windows compatibilty
%               try/catch renaming
%               Use java as it is blind to your system.
%               Robust path additions
%               Preference adding.
%               Logging question
%               New getpref/setpref structure

%% Get system information
% Use Java, it's more reliable!
st_home = cast(java.lang.System.getProperty('user.home'),'char');
st_user = cast(java.lang.System.getProperty('user.name'),'char');
try
    st_comp = cast(java.net.InetAddress.getLocalHost.getHostName,'char');
catch
    st_comp = 'unavailable';
end

%% Add libraries to the path
% ! NOTE ! We do not have to add files to path if this is correct!
% Go to default place ???
try % Try the new way
    libroot = sdext.getpref('libroot').val;
catch % We have not fixed the path
    try
        libroot = getpref('mtools','libroot'); % Try the old way
        addpath(genpath(fullfile(libroot,'Spectra')))
    catch
        if ismac
            libroot = uigetdir(fullfile(matlabroot,'toolbox'),'Select mtools root directory');
        else
            libroot = uigetdir(st_home,'Select mtools root directory');
        end
        if all(libroot == 0)
            error('Without an mtools directory, here be dragons!')
        else
            addpath(genpath(fullfile(libroot,'Spectra')))
            sdext.setpref('libroot',libroot)
            choice = questdlg('Would you to enable experimental features?', ...
                'Enable Extras', ...
                'Yes','No','No');
            % Handle response
            switch choice
                case 'No'
                    sdext.setpref('experimental',0)
                case 'Yes'
                    sdext.setpref('experimental',1)
            end
        end
    end
end

%% Start a diary
doLog = sdext.getpref('doLog').val;

if (doLog) == 2
    choice = questdlg('Would you to enable logging?', ...
        'Enable logging', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'No'
            sdext.setpref('doLog',0)
            doLog = 0;
        case 'Yes'
            sdext.setpref('doLog',1)
            doLog = 1;
    end
end

d = fullfile(st_home,'Documents','MATLAB');
if isdir(d)
    cd(d)
else
    d = [getenv('HOME') filesep 'MATLAB'];
    if isdir(d)
        cd(d)
    end
end

if doLog
    if exist(fullfile(d,'matlab.log'),'file') == 2
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
    diary(fullfile(d,'matlab.log'))
    disp([ 'matlab.log started on : ' datestr(now) ])
end

%% User defined startup commands
% This file contains user commands. i.e setting up other programs
if exist('startuser.m','file')
    eval('startuser');
end

if all(~[isempty(st_user)  isempty(st_comp)])
    fprintf( 'Welcome %s@%s! Spectra files are loaded\n',st_user,st_comp);
end

%% Fix the broken stuff!
% If you have the curve fitting toolbox horace and herbert overwrite some functions
v = ver;
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

%% Load previous workspace
if exist('matlab.mat','file')
    disp('Getting previous workspace data from matlab.mat')
    load matlab.mat
end