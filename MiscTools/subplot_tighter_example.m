%% subplot_tight example script
% This is a brief tutorial of subplot_tighter. For better control of subplots
% Note: You can not use the 'align' argument as it is already aligned!
%       If your lables are outside of the plot you have to change Vaules 1-4.

%% Make a referece
% Value 1 = Top and bottom margins
% Value 2 = Left and right margins
% Value 3 = Vertical inter-plot spacing
% Value 4 = Horizontal inter-plot spacing
margins = [0.01 0.01 0.2 0.1];

%% Example 1
% A 2 by 2 grid with no spread
figure
subplot_tighter(2,2,1,margins)
subplot_tighter(2,2,2,margins)
subplot_tighter(2,2,3,margins)
subplot_tighter(2,2,4,margins)

%% Example 2
% A 2 by 2 grid with vertical spread
figure
subplot_tighter(2,2,[1 3],margins)
subplot_tighter(2,2,2,margins)
subplot_tighter(2,2,4,margins)

%% Example 3
% A 2 by 2 grid with horizontal spread
figure
subplot_tighter(2,2,[1 2],margins)
subplot_tighter(2,2,3,margins)
subplot_tighter(2,2,4,margins)

%% Example 4
% In this example we will have 2 sets of margins. 
% !! Note that Value 1 and Value 2 have to be the same for both margins !!
figure
margins = [0.05 0.05 0.075 0.05];
subplot_tighter(3,3,1,margins)
subplot_tighter(3,3,2,margins)
subplot_tighter(3,3,3,margins)
margins = [0.05 0.05 0.0 0.05];
subplot_tighter(3,3,4,margins)
subplot_tighter(3,3,5,margins)
subplot_tighter(3,3,6,margins)
subplot_tighter(3,3,7,margins)
subplot_tighter(3,3,8,margins)
subplot_tighter(3,3,9,margins)

