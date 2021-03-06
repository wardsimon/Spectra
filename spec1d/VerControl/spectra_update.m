function onlineRev = spectra_update(installDir)
% updates the Spectra installation from the internet
%
% SW_UPDATE()
%
% {onlineRev} = SPECTRA_UPDATE(installDir)
%
% spectra_update creates a new folder with the latest release beside the current
% Spectra installation and add the new version to the working path (and
% removing the old one).
%
% Input:
%
% installDir    Folder name, where the new version is installed. Default is
%               the parent folder of the current version of Spectra. If
%               installDir == '.' update will be installed to current
%               folder.
%
% Output:
%
% onlineVer     If output is expected, the revision number of the online
%               Spectra is given. Optional.
%
% This has been modified from the spinW help file. Thanks to Sandor Toth!

% check current version
swVer = spectra_version;

% base url, where the sw_download_info file stored
baseUrl = 'https://raw.githubusercontent.com/simonward86/Spectra/LNS/spec1d/VerControl/this_version.txt';

if nargout == 0
    if ~isempty(swVer.Version)
        answer = getinput(...
            ['This is a not yet released version of Spectra, update is not recommended!\n'...
            'Do you want to continue? (y/n)'],'yn');
        
        if answer == 'n'
            disp('Spectra update process cancelled!');
            return
        end
        swVer.Revision = 0;
    end
    
    if nargin == 0
        installDir = spectra_rootdir;
        strIdx = strfind(installDir,filesep);
        installDir = installDir(1:strIdx(end-1));
    else
        if installDir(1) == '.'
            installDir = [pwd filesep];
        elseif installDir(end) ~= filesep
            installDir = [installDir filesep];
        end
    end
end

% download the link to the newest version & comments!
% the file format:
%   link to new download
%   revision number
%   message in the next few lines
try
    newInfo = urlread(baseUrl);
catch
    error('spectra_update:NoNetwork','It looks like there is a problem with your network connection!');
end

% separate lines of text
newInfo = textscan(newInfo, '%s', 'delimiter', sprintf('\n'));
newInfo = newInfo{1};

newLink = newInfo{1};
newRev  = str2double(newInfo{2});

if nargout == 1
    % Give the revision number and exit.
    onlineRev = newRev;
    return
end

% check whether the online version is newer (compare revision numbers)
if ischar(swVer.Revision)
    swVer.Revision = str2double(swVer.Revision);
end

if swVer.Revision ==  newRev
    disp('Your Spectra installation is up to date!');
    return
elseif swVer.Revision >  newRev
    disp('Your Spectra installation is newer than the one in the repository! Lucky you! :)');
    return
end

fprintf('Current version has a revision number: %d\n',swVer.Revision);
fprintf('New version has a revision number:     %d\n',newRev);

answer = getinput('Do you want to continue? (y/n)','yn');

if answer == 'n'
    disp('Spectra update process cancelled!');
    return
end

fprintf('New version will be installed to: %s\n',installDir);
answer = getinput('Do you want to continue? (y/n)','yn');

if answer(1) == 'n'
    disp('Spectra update process cancelled!');
    return
end

% save new update as a zip file
updateName = 'Spectra_update_files.zip';
fprintf('Downloading update from %s... ',newLink);
urlwrite(newLink,[installDir updateName]);
fprintf('ready!\n');

% decompress zip file
zipList = unzip([installDir updateName],installDir);

% get folder name
folName = [installDir strtok(zipList{1}(numel(installDir)+1:end),filesep)];

% remove old Spectra installation from path
fprintf('\nRemoving path to old Spectra installation!\n')
rmpath(genpath(spectra_rootdir));

% adding new path
fprintf('Adding path to new Spectra installation: %s!\n',folName);
ww = warning;
warning('off');
addpath(genpath(folName));
warning(ww);

answer = getinput('Do you want to save the new path (savepath)? (y/n)','yn');

if answer == 'y'
    savepath
    disp('The current path is saved!');
end

disp('Removing unnecessary files... ')
delete([installDir updateName]);

% Below section is implemented in install_spinw
% % remove files aren't needed for new Matlab versions
% % functions introduced in R2014a
% if ~verLessThan('matlab', '8.1')
%     % strjoin()
%     fList = dir([folName filesep 'external' filesep 'strjoin*']);
%     for ii = 1:numel(fList)
%         delete([folName filesep 'external' filesep fList(ii).name]);
%     end
%     % strsplit
%     fList = dir([folName filesep 'external' filesep 'strsplit*']);
%     for ii = 1:numel(fList)
%         delete([folName filesep 'external' filesep fList(ii).name]);
%     end
%     % gobjects
%     fList = dir([folName filesep 'external' filesep 'gobjects*']);
%     for ii = 1:numel(fList)
%         delete([folName filesep 'external' filesep fList(ii).name]);
%     end
% else
%     % rename the functions to be used in Matlab versions prior to R2014a
%     % strjoin()
%     try %#ok<*TRYNC>
%         movefile([folName filesep 'external' filesep 'strjoin0.m'],...
%             [folName filesep 'external' filesep 'strjoin.m']);
%     end
%     % strsplit()
%     try
%         movefile([folName filesep 'external' filesep 'strsplit0.m'],...
%             [folName filesep 'external' filesep 'strsplit.m']);
%     end
%     % gobjects()
%     try
%         movefile([folName filesep 'external' filesep 'gobjects0.m'],...
%             [folName filesep 'external' filesep 'gobjects.m']);
%     end
% end
% 
% % functions introduced in R2015a
% % if ~verLessThan('matlab', '8.5')
% %     % uniquetol()
% %     fList = dir([folName filesep 'external' filesep 'uniquetol*']);
% %     for ii = 1:numel(fList)
% %         delete([folName filesep 'external' filesep fList(ii).name]);
% %     end
% % end
% 
% fprintf(['In order to reach SpinW after restarting Matlab, the following\n'...
%     'line has to be added to your startup.m file:\n']);
% fprintf('  addpath(genpath(''%s''));\n',folName);
% 
% % location of Matlab startup file
% sfLoc = which('startup');
% uPath = userpath;
% % remove ':' and ';' characters from the userpath
% uPath = [uPath(~ismember(uPath,':;')) filesep 'startup.m'];
% 
% % create new startup.m file
% if isempty(sfLoc)
%     answer = getinput(sprintf(['You don''t have a Matlab startup.m file,\n'...
%         'do you want it to be created at %s? (y/n)'],uPath),'yn');
%     if answer == 'y'
%         fclose(fopen(uPath,'w'));
%         sfLoc = uPath;
%     end
% end
% 
% if ~isempty(sfLoc)
%     
%     answer = getinput(sprintf(['Would you like to add the following line:\n'...
%         sprintf('addpath(genpath(''%s''));',folName) '\nto the end of '...
%         'your Matlab startup file (%s)? (y/n)'],sfLoc),'yn');
%     
%     if answer == 'y'
%         fid = fopen(sfLoc,'a');
%         fprintf(fid,['\n%%###SW_UPDATE\n%% Path to the SpinW (rev. %d) '...
%             'toolbox:\naddpath(genpath(''%s''));\n%%###SW_UPDATE\n'],newRev,folName);
%         fclose(fid);
%     end
% end

%% This is for code options. Not implemented yet
% if numel(newMsg)>0
%     answer = getinput('Do you want to see the release information? (y/n)','yn');
%     if answer == 'y'
%         disp('Release information:')
%         disp(repmat('-',[1 60]))
%         for ii = 1:numel(newMsg)
%             fprintf('\t%s\n',newMsg{ii});
%         end
%         disp(repmat('-',[1 60]))
%     end
% end

% answer = getinput(...
%     ['\nIn order to refresh the internal class definitions of Matlab (to\n'...
%     'access the new SpinW version), issuing the "clear classes" command\n'...
%     'is necessary. However this command will also clear all your variables\n'...
%     'in the Matlab internal memory. Would you like the updater to issue\n'...
%     'the command now, otherwise you can do it manually later.\n'...
%     'Do you want to issue the command "clear classes" now? (y/n)'],'yn');
% 
% if answer == 'y'
%     clear('classes'); %#ok<CLCLS>
%     disp('Matlab class memory is refreshed!')
% end

disp('Update was successful!')

answer = getinput('Do you want to run the spectra_on command from the update? (y/n)','yn');

switch answer
    case 'n'
        disp('Don''t forget to run spectra_on later!');
    case 'y'
        spectra_on;
end

end

function answer = getinput(message,good)
% get the necessary letter input

answer = ' ';
while ~ismember(answer(1),good)
    answer = input(message,'s');
    if isempty(answer)
        answer = 0;
    end
end
answer = answer(1);

end