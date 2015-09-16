function varargout=ALPS_H(varargin)

global g

if isempty(g)
    g=2.06;
end

if nargin>2
    for i=1:2:(nargin)-1
        filenames{i}=varargin{i};
        ptfi=1:varargin{i+1};
    end
    ldt=varargin{nargin};
else
    filenames=varargin{1};
    ldt='magnetization';
    ptfi=1;
end

r=1;
for j=1:length(filenames)
    filename=filenames{j};
    for i=ptfi
        temp=xmlread([filename '.task' num2str(i) '.plot.' ldt '.xml']);
        temp2=temp.getChildNodes;
        temp3=temp2.item(1);
        allListitems=temp3.getElementsByTagName('point');
        M=zeros(allListitems.getLength,2);
        for k = 0:allListitems.getLength-1
            thisListitem = allListitems.item(k);
            thisList = thisListitem.getElementsByTagName('y');
            thisElement = thisList.item(0);
            M(k+1,2) = str2double(char(thisElement.getFirstChild.getData));
            pn=thisListitem.getParentNode;
            temp = char(pn.getAttribute('label'));
            M(k+1,1)=1.4887*str2double(temp(3:end))/g;
        end
        varargout{r}=M;
        r=r+1;
    end
end