function varargout=resortALPS(varargin)

k=1;
name_given=0;
for i=1:nargin
    if isa(varargin{i},'ALPS')
        for j=1:length(varargin{i})
            ALPS(k)=varargin{i}(j);
            k=k+1;
        end
    else
        if isa(varargin{i},'cell')
            x_name=varargin{i}{1};
            y_name=varargin{i}{2};
            name_given=1;
        end
    end
end

if ~name_given
    error('Please supply your x and y names')
end

for i=1:(k-1)
    if strcmp(x_name,ALPS(i).xname)
        if strcmp(y_name,ALPS(i).zname)
            temp=ALPS(i).y;
            tempN=ALPS(i).yname;
            ALPS(i).y=ALPS(i).z;
            ALPS(i).z=temp;
            ALPS(i).yname=ALPS(i).zname;
            ALPS(i).zname=tempN;
        end
    elseif strcmp(x_name,ALPS(i).yname)
        temp=ALPS(i).x;
        tempN=ALPS(i).xname;
        ALPS(i).x=ALPS(i).y;
        ALPS(i).xname=ALPS(i).yname;
        if strcmp(y_name,ALPS(i).xname)
            ALPS(i).y=temp;
            ALPS(i).yname=tempN;
        elseif strcmp(y_name,ALPS(i).zname)
            ALPS(i).y=ALPS(i).z;
            ALPS(i).z=temp;
            ALPS(i).yname=ALPS(i).zname;
            ALPS(i).zname=tempN;
        else
            warning('No uncomplete sorting applied, check names')
        end
    elseif strcmp(x_name,ALPS(i).zname)
        temp=ALPS(i).x;
        ALPS(i).x=ALPS(i).z;
        tempN=ALPS(i).zname;
        ALPS(i).xname=tempN;
        if strcmp(y_name,ALPS(i).xname)
            ALPS(i).z=ALPS(i).y;
            ALPS(i).zname=ALPS(i).yname;
            ALPS(i).y=temp;
            ALPS(i).yname=tempN;
        elseif strcmp(y_name,ALPS(i).yname)
            ALPS(i).z=temp;
            ALPS(i).zname=tempN;
        else
            warning('No uncomplete sorting applied, check names')
        end
    else
        warning('No sorting applied, check names')
    end
end

varargout{1}=ALPS;