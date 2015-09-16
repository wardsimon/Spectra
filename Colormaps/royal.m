function cmap = royal(n)
%   colormap(royal)
%
%   Initializes custom colourmap. Colours are: black, 
%   Royal blue (hence the name), cerulean, turquois, 
%   green, yellow, white.
%
%   Composed by Matt Mena, 17.07.2014


% Set the main colours.
mattmapp=zeros(5,3);
mattmapp(2,:)=[0 0.137 0.340]; % Royal blue
mattmapp(3,:)=[0 0.481 0.652]; % Cerulean
mattmapp(4,:)=[0.25 0.875 0.8125]; % Turquois
mattmapp(5,:)=[0.5 0.9 0.4]; % Green
mattmapp(6,:)=[1 1 0]; % Yellow
mattmapp(7,:)=[1 1 1]; % white

nmappoi=20;

% Default size is [121,3], as all other colormaps
cmap=zeros(nmappoi*(numel(mattmapp(:,1))-1)+1,3);

for j=1:numel(mattmapp(:,1))-1
    mapploc=zeros(20,3);
    for k=1:nmappoi
        mapploc(k,:)=mattmapp(j,:)*(nmappoi-k+1)/nmappoi+...
            mattmapp(j+1,:)*(k-1)/nmappoi;
    end
    cmap(((j-1)*nmappoi+1):(j*nmappoi),:)=mapploc;
end

cmap(end,:)=mattmapp(end,:);

% Assuming that someone asked for a specific number of color levels, here
% they get computed
if nargin==1
    % Create temporary cmap2
    cmap2 = zeros(n,3);
    cmap2(1,:) = cmap(1,:);
    for j=2:n-1
        indl = (j-1)*120/(n-1)+1;
        dk = mod(indl,1); k = indl-dk;
        col = cmap(k,:)*(1-dk)+cmap(k+1,:)*dk;
        cmap2(j,:) = col;
    end
    cmap2(end,:) = cmap(end,:);
    % Substitute cmap2 to cmap.
    cmap = cmap2;
end
