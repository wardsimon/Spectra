function update_horace(varargin)

    try
        horace_off
        hor_run = 1;
        fprintf('Turning Horace off\n')
    catch
        fprintf('Horace is off\n')
        hor_run = 0;
    end
   
    libroot = sdext.getpref('libroot').val;
    
    if nargin == 0
        url = 'http://svn.isis.rl.ac.uk/viewvc/Horace/trunk/?view=tar';
    else
        if ischar(varargin{1})
            url = varargin{1};
        end
    end
    
    if exist(fullfile(libroot,'Horace'),'file')
        [fstatus, message] = copyfile(fullfile(libroot,'Horace'),fullfile(libroot,'Horace_Backup'));
    else
        fstatus = 2;
    end
    
    if fstatus == 0
        error('Failed to backup Horace \n%s',message)
    end
    
    fprintf('Getting update from:\n%s\nThis may take some time....\n',url)
    [filestr,status] = urlwrite(url,fullfile(libroot,'update.tar.gz'));
    
    if status~=1
        error('Failed to download Horace \n%s',message)
    else
        fprintf('Download succeeded.\n')
    end
    
    if fstatus ~= 2
        [status, message]= rmdir(fullfile(libroot,'Horace'),'s');
    else
        status = 1;
    end
    
    if status~=1
        error('Failed to remove previous Horace version\n%s',message)
    end
    
    files = untar(filestr);
    movefile(fullfile(libroot,fileparts(files{1}),'*'),fullfile(libroot,'Horace'))
    
    rmdir(fullfile(libroot,fileparts(files{1})))
    
    delete(filestr)
    
    fprintf('Sucessfully updated Horace\n')
    
    if hor_run
        horace_on
    end