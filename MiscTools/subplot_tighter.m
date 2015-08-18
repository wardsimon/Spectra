function h=subplot_tighter(m,n,p,margins,varargin)
%function subplot_tighter(m,n,p,margins,varargin)
%
% Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the margins between
% neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the margins between
% subplots can reach 40% of figure area, which is pretty lavish.  
%
% Input arguments (defaults exist):
%   margins- four elements vector [vertical,horizontal vertical_interplot horizontal_interplot] defining the margins between neighbouring axes. Default value
%            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
%            relatively large axis. 
%
% Output arguments: same as subplot- none, or axes handle according to function call.
%
% Issues & Comments: Note that if additional elements are used in order to be passed to subplot, margins parameter must
%       be defined. For default margins value use empty element- [].      
%
% Original author and Date:  Nikolay S. 29/03/2011. 
% Last update             :  Ward S. 13/05/2014

if (nargin<4) || isempty(margins)
    margins=[0.01,0.01, 0.01, 0.01]; % default margins value- 1% of figure
end

if length(margins)==1
    margins(2:4)=margins;
end

if ~iscell(margins)
    margins = {{margins(1) margins(1)},{margins(2) margins(2)},margins(3), margins(4)};
end

%note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
[subplot_col,subplot_row]=ind2sub([n,m],p);  


height = (1 - (margins{1}{1} + margins{1}{2}) -(m-1)*margins{3})/m;  % single subplot height
width  = (1 - (margins{2}{1} + margins{2}{2}) -(n-1)*margins{4})/n;  % single subplot width

% note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
subplot_cols = 1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
subplot_rows = 1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

merged_height = subplot_rows*( height+margins{3} )- margins{3};   % merged subplot height
merged_width  = subplot_cols*( width +margins{4} )- margins{4};   % merged subplot width

merged_bottom = (m-max(subplot_row))*(height+margins{3}) + margins{1}{1}; % merged subplot bottom position
merged_left   = margins{2}{1} + (min(subplot_col)-1)*(width + margins{4}); % merged subplot left position

pos_vec = [merged_left merged_bottom merged_width merged_height];

% h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
% Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
h_subplot = subplot('Position',pos_vec,varargin{:});

if nargout~=0
    h=h_subplot;
end