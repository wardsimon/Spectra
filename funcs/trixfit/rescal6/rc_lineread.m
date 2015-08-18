function [x]=lineread(file,par)
    %
    % function [x]=lineread(file,par)
    % Read parameters from ILL files.
    x=[];
    err=[];
    
    %------------------ Open data file --------------------------------
    com=[1:400];
    fid=fopen(file,'r');
    if (fid<0)
        error(['file ' file ' not found']);
    end
    
    while ~feof(fid)
        com = fgetl(fid);
        if (findstr(par,com)==1)
            streq    = findstr('=',com);
            strcomma = findstr(',',com);
            for i=1:(max(size(streq))-1)
                x = [x str2double(com(streq(i)+1:strcomma(i)-1))];
            end
            if (max(size(streq))>1)
                x = [x str2num(com(streq(i+1)+1:max(size(com))-1))];
            end
            if isempty(strcomma)
                x = [x str2double(com((streq+2):end))];
            end
%                 x(isnan(x))=[];
        end
    end
    
    fid=fclose(fid);
    return
