function out=plus(varargin)

k=1;
for i=1:nargin
    if isa(varargin{i},'ALPS')
        for j=1:length(varargin{i})
            if i==1 & k==1
                in=varargin{i}(j);
            else
                ALPS(k)=varargin{i}(j);
                k=k+1;
            end
        end
    end
end

if k>1
    x=in.x;
    y=in.y;
    z=in.z;
else
    error('Plese supply more than 1 ALPS object')
end

for i=1:length(ALPS)
    xnew=ALPS(i).x;
    ynew=ALPS(i).y;
    znew=ALPS(i).z;
    
    out(i)=in;
    
    if sum(x==xnew)==length(x)
        out(i).y=y+ynew;
        out(i).z=z+znew;
    else
        out(i).y=y+interp1(xnew,ynew,x);
        out(i).z=z+interp1(xnew,znew,x);
    end
end

% out=class(out,'ALPS');