function [xx, yy, ee]=mv_rebin(x, y, err, dx)
%
% function [xx, yy, ee]=mv_rebin(x, y, err, dx)
%
% Rebin x, y, err data in bins specified by dx.
% If dx is a vector, it specifies the bin boundaries,
% otherwise it specifies the bin widths.
%
%  MZ 28.3.95

% Construct bin boundaries if necessary
if length(dx)==1        
	dx=min(x):dx:max(x);
end

% Construct new x-vector (bin centres)
xx=mean([dx(1:length(dx)-1); dx(2:length(dx))]);
yy=zeros(size(xx));
ee=yy;

% Fill each bin in turn
for i=1:length(xx)
	j=find((x>=dx(i)) & (x<dx(i+1)));
	yy(i)=sum(y(j))/length(j);
	ee(i)=sqrt(sum(err(j).^2))/length(j);
end
xx=xx';
yy=yy';
ee=ee';
