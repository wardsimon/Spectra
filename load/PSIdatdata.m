function [x, y, e, xname, yname,mname,optpars]=PSIdatdata(filespec)
    % Load .dat data from PSI files.
    % Works with:
    %   Trics
    %   Morpheus
    %   Orion
    %   Eiger
    %   Tasp
    % Simon Ward May 2011, simon.ward@psi.ch
    
    x=[]; y= []; e=[]; yname=''; xname=''; mname = '';
    
    %----- Parse filespec --------------------------------------
    
    [fspec, filespec] = strtok(filespec,',');
    
    while ~isempty(filespec)
        [s, filespec] = strtok(filespec,',');
        fspec = str2mat(fspec,s);
    end
    
    [nargs, nchars] = size(fspec);
    
    %----- Update scan parameters from filespec---------------------------
    
    xname = nameit('X=');
    yname = nameit('Y=');
    mname = nameit('M=');
    
    filename = deblank(fspec(1,:));
    
    
    % ---- Get Data from file ----------------------------------
    fid = fopen(filename,'r');
    finishup = onCleanup(@() fclose(fid));
    
    if fid == -1
        warning('Matlab:FileNotFound','The file is not found!')
        return
    end
    
    i = 1;
    
    while ~feof(fid)
        linDat = fgetl(fid);
        if strmatch('NP',linDat) | strmatch('PNT',linDat)
            [lspec, linspec] = strtok(linDat,' ');
            while ~isempty(linspec)
                [s, linspec] = strtok(linspec,' ');
                lspec = str2mat(lspec,s);
            end
            for j = 1:length(lspec(:,1))
                if strcmp(deblank(lspec(j,:)),xname)
                    xpos = j;
                end
                if strcmp(deblank(lspec(j,:)),yname)
                    ypos = j;
                end
                if strcmp(deblank(lspec(j,:)),mname)
                    mpos = j;
                end
            end
            
            tdata = importdata(filename,' ',i+1);
            data  = tdata.data;
            x = data(:,xpos);
            y = data(:,ypos)./data(:,mpos);
            e = sqrt(data(:,ypos))./data(:,mpos);
            break;
        end
        i = i+1;
    end
    
    function name  = nameit(f)
        k = strmatch(f,fspec);
        if ~isempty(k)
            name = deblank(fspec(k(end),3:nchars));
        else
            error('PSIdatdata:VairableNotFound','The vairable %s is not found.',f)
        end
    end
end