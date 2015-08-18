function [Harrow]=gf_arrow(X,Y,Arrcol,lwidth,Arrtxt,txtoffset)
% MATLAB graphics function to create an arrow object
%
%
% function [Harrow]=gf_arrow(X,Y,Arrcol,lwidth,Arrtxt)
% [Harrow] is the array of handles [Line; Flight; Text]
% [X] is the start and end x-value.
% [Y] is the start and finish y-value.
% [Arrcol] is the RGB colour value of arrow and label {default white}
% lwidth is the line width {default 1}
% [Arrtxt] label that is put on the rhs of the arrow
%
% DAT 21.11.97 ORNL

if (nargin<6) txtoffset=0.1; end
if (nargin<5) Arrtxt=['']; end
if (nargin<4) lwidth=1; end 
if (nargin<3) Arrcol=[1 1 1]; end


% draw line
Hline=line(X,Y);
% draw flights on the arrow indicating direction
% this is composed of a triple of points
% We shall make the size of the flight proportianate
% with the arrow length
Darrow=[diff(X) diff(Y)];
Carrow=[sum(X)/2 sum(Y)/2];
prop=.05;
% generate three points which are at a distance prop*Darrow
% away from Carrow at angles of {0,1,2}*pi/3
theta=2*pi/3;
% use a rotation operator to generate the points
R=[cos(theta) sin(theta); -sin(theta) cos(theta)];
Fpt=[];
for i=-1:1
  Fpt=[Fpt ; Carrow+prop*(Darrow*(R^i)')];
end
Hflight=line(Fpt(:,1),Fpt(:,2));
txtpos=Carrow+2*txtoffset*Darrow*[0 1; -1 0]';
Htxt=text(txtpos(1),txtpos(2),Arrtxt);

Harrow=[Hline;Hflight;Htxt];
set(Harrow(1:2),'color',Arrcol,'linewidth',lwidth)
set(Harrow(3),'color',Arrcol);

return
end
% -------------- end of arrow formation ----------------------
