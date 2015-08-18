function of=ALPS(varargin)

r=1;
for i=1:nargin
    for j=1:length(varargin{i}.measurements)
        for k=1:varargin{i}.tasks
            % Get to the data store
            temp=xmlread([varargin{i}.filename '.task' num2str(k) '.plot.' varargin{i}.measurements{j} '.xml']);
            temp2=temp.getChildNodes;
            temp3=temp2.item(1);
            
            % Read some data
            of(r).observable=char(temp3.getElementsByTagName('PARAMETER').item(0).getFirstChild.getData);
            of(r).xname=char(temp3.getElementsByTagName('xaxis').item(0).getAttribute('label'));    of(r).xname=of(r).xname(1);
            of(r).yname=char(temp3.getElementsByTagName('yaxis').item(0).getAttribute('label'));    of(r).yname=of(r).yname(1);
            tl=char(temp3.getElementsByTagName('set').item(0).getAttribute('label'));
            n=strfind(tl,'=');
            of(r).zname=tl(1:(n-1));
            
            % Extract the simulation data
            x=[]; y=[];z=[];
            for l=0:(temp3.getElementsByTagName('set').getLength-1);
                tl=char(temp3.getElementsByTagName('set').item(l).getAttribute('label'));
                for m = 0:(temp3.getElementsByTagName('set').item(l).getElementsByTagName('point').getLength-1)
                    x(end+1)=str2double(char(temp3.getElementsByTagName('set').item(l).getElementsByTagName('point').item(m).item(0).getFirstChild.getData));
                    y(end+1)=str2double(char(temp3.getElementsByTagName('set').item(l).getElementsByTagName('point').item(m).item(1).getFirstChild.getData));
                    z(end+1)=str2double(tl((n+1):end));
                end
            end
            
            % Correct field from K into T
            if strcmp(of(r).xname,'h')
                of(r).x=1.4887*x/varargin{i}.g;
            else
                of(r).x=x;
            end
            if strcmp(of(r).yname,'h')
                of(r).y=1.4887*y/varargin{i}.g;
            else
                of(r).y=y;
            end
            if strcmp(of(r).zname,'h')
                of(r).z=1.4887*z/varargin{i}.g;
            else
                of(r).z=z;
            end
            r=r+1;
        end
    end
end

of=class(of,'ALPS');