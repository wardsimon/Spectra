function varargout = colorwheel(varargin)
    % Create a colorwheel in 1-3D.
    %
    % The options are:
    %     options.steps = 15; -> Number of points in unknown dimension
    %     options.h = [];     -> If you want to specify H
    %     options.s = [];     -> If you want to specify S
    %     options.v = [];     -> If you want to specify V
    %     options.plot = 0;   -> Do you want to plot the colorwheel
    %     options.sectors = 4;-> Often you don't want colors next to each other. This is how much spread.
    %     options.phase_offset = 240;-> Try to get blue as #1!
    %
    % The output is a [options.steps 3] matrix of rgb values sorted by options.sectors if applicable.
    %
    % Simon Ward 24/02/2014
    % simon.ward@psi.ch
    
    
    if nargin == 0 
        options = struct();
    else
        options = varargin;
    end
    
    %parse inputs
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    p.addOptional('h',[],@(x) all(x<=1 & x>=0));
    p.addOptional('s',[],@(x) all(x<=1 & x>=0));
    p.addOptional('v',[],@(x) all(x<=1 & x>=0));
    p.addOptional('plot',0,@isnumeric);
    p.addOptional('sectors',4,@isnumeric);
    p.addOptional('phase_offset',240,@isnumeric);
    p.addOptional('steps',15,@isnumeric);
  
    if iscell(options)
        p.parse(options{:});
    elseif isstruct(options)
        p.parse(options);
    else
        error('Can not parse input')
    end
    options = p.Results;
 
    options.phase_offset = deg2rad(options.phase_offset);
    
    base = linspace(0,1,options.steps+1);
    j = 0;
    io = [0 0 0];
    if isempty(options.h)
        j = j+1;
        io(1) = 1;
    end
    if isempty(options.s)
        j = j+1;
        io(2) = 1;
    end
    if isempty(options.v)
        j = j+1;
        io(3) = 1;
    end
    
    % Generate matrices
    C = cell(j,1);
    [C{:}] = ndgrid(base);
    doit = cellfun(@(x){x(:)}, C);
    doit = [doit{:}];
    j = 1;
    hsv_matrix = zeros(length(doit(:,1)),3);
    if io(1) == 0
        hsv_matrix(:,1) = options.h;
    else
        hsv_matrix(:,1) = doit(:,j);
        j = j + 1;
    end
    if io(2) == 0
        hsv_matrix(:,2) = options.s;
    else
        hsv_matrix(:,2) = doit(:,j);
        j = j + 1;
    end
    if io(3) == 0
        hsv_matrix(:,3) = options.v;
    else
        hsv_matrix(:,3) = doit(:,j);
    end
    
    
    rgb_matrix = hsv2rgb(hsv_matrix);
    rgb_matrix = rgb_matrix(1:(end-1),:);
    
    theta = 2*pi*hsv_matrix(1:(end-1),1);
    
    % Plot the colorwheel if needed
    if options.plot
        
        r = hsv_matrix(1:(end-1),2);
        [x, y] = pol2cart(theta,r);
        z = hsv_matrix(1:(end-1),3);
        
        figure;
        hold on
        for i = (1:(numel(x)))
            plot3(x(i), y(i),z(i),'MarkerFaceColor',rgb_matrix(i,:),'Marker','o','MarkerEdgeColor','None');
        end
        [x y] = pol2cart(deg2rad(0:0.1:360),1.05);
        plot(x,y,'k-')
        offset = 0;
        if options.plot > 1
            offset = options.phase_offset;
        end
        for i = 1:options.sectors
            [x, y] = pol2cart(((i-1)*(2*pi/options.sectors)) + offset,1.05);
            plot([0 x],[0 y],'k-')
            [x, y] = pol2cart(mean(((i-(1:2))*(2*pi/options.sectors)) + offset),1.05);
            text(x,y,sprintf('Sector %i',i))
        end
        axis([-1.1 1.1 -1.1 1.1])
    end
    
    % Apply a phase offset (i.e, blue first :-) )
    tmp = abs(theta-options.phase_offset);
    [~, ind] = min(tmp);
    rgb_matrix = circshift(rgb_matrix,[ind-2 0]);
    theta = circshift(theta,[ind-2 0]);
    
    % Re-sort the colors if need be
    if length(theta) > 1
        if diff(theta(1:2)) < 2*pi/options.sectors
            ind = mod(0:(length(theta)-1),options.sectors);
            rgb_matrix2 = [];
            for i = 0:(options.sectors-1)
                rgb_matrix2 = vertcat(rgb_matrix2,rgb_matrix(ind == i,:));
            end
            rgb_matrix = rgb_matrix2;
        end
    end
    % Create outputs
    if nargout == 1
        varargout{1} = rgb_matrix;
    end
end
function rad = deg2rad(deg)

% deg2rad (rad)
%
% Keith Brady 18-08-2010
%
% Conversion function from Degrees to Radian.
% an alterative form to the functions by Richard Medlock renamed as
% suggested in the blog article 
%
% http://blogs.mathworks.com/pick/2009/11/27/degrees-and-radians/

rad = (pi/180).*deg;
end
function deg = rad2deg(rad)

% rad2deg (rad)
%
% Keith Brady 18-08-2010
%
% Conversion function from Radians to Degrees.
% an alterative form to the functions by Richard Medlock renamed as
% suggested in the blog article 
%
% http://blogs.mathworks.com/pick/2009/11/27/degrees-and-radians/


deg = rad./(pi/180);
end