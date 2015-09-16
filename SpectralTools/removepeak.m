function [nx,ny,iy,sy] = removepeak(xs,ys, thr, pos,width)
% removepeak : Removes some peaks in a signal
%Syntax: [new_x,new_y,interpolated_y,smoothed_y] = removepeak(x,y,{thr=1, pos, width}) or new_y = removepeak(y)
%
% Removes some peaks in 'y', with user prompt.
% Returns a shorter signal for x and y axis (peaks data removed), the original
% signal with linear interpolation instead of removed peaks, and smoothed y.
% optional : thr threshold for peaks searching. Usually from 0.1 to 10.
%            pos, width = part of 'y' to be scanned for peaks around pos.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  removes some peaks in a signal

% uses : findpeaks.m
% Spectral Tools, E.Farhi 07/96

if (nargin < 1)
	error ('usage : new_y = remove_peak(x,y,{thr, pos, width}) or new_y = remove_peak(y)');
end

if (nargin == 1)
	ys=xs;
	xs=1:length(ys);
end

x=xs(:);
y=ys(:);

nd = length(x);

if ~exist('width')
	width = round(max(3,nd / 50));
end

if isempty(width)
	width = round(max(3,nd / 50));
end

if ~exist('pos') | ~exist('width')
	pos = x(ceil(nd/2));
	width = abs(x(nd)-x(1));
end

if isempty(pos) | isempty(width)
	pos = x(ceil(nd/2));
	width = abs(x(nd)-x(1));
end

if ~exist('thr')
	thr = 1;
end

if isempty(thr)
	thr = 1;
end

t=find( (x>=(pos-width)) & (x<=(pos+width)) );		% index zone around pos
disp('Looking for peaks...');
y = y(t);
x = x(t);
iy = y;
wt = ones(1,length(y));

flag = 0;
while (flag == 0)		% choosing peak to remove, and confirm choice.

  [peakmat,pmax,sy,ny] = findpeaks(x,y,thr,[],[],'silent');
  fprintf(1,'Found %i peaks.\n', rows(peakmat));
  fprintf(1,'Noise is about %f .\n',ny);
  clg;
  title('Peaks position in signal');
  ylabel('Y');
  xlabel('X');
  plot(x,y,x,pmax(:).*y)
  if (~isempty(peakmat))
	fprintf(1,' Number  Position  Intensity  (ExpValue)  Width      Index\n') ;
	for index = 1:rows(peakmat)
%          [ index left_width right_width  max_pos Intensity Width ].
	   fprintf(1,'Peak %3i X=%7.2f I=%8.2g (%8.2g) W1=%7.2f  i=%i\n', index, peakmat(index,4), peakmat(index,5), y(peakmat(index,1)), peakmat(index,6), peakmat(index,1));
	end
	flag = 1;
  else
	thr = thr / 2;
  end
end
flag = 0;
while (flag == 0)
	disp(' ');
	index = input('Enter peak number to remove (or Enter to quit) : ');
	if (isempty(index))
		flag = 2;	% quit
	else 
		flag = 0;
				% replace peak by noise value and interpolation
		S=peakmat(index,4);
		W=peakmat(index,6);
		id = peakmat(index,1);
		WL = mean([ abs(x(max(1,id-peakmat(index,2)))-x(id)) W/1.4 ]);
		WR = mean([ abs(x(min(length(x),id+peakmat(index,3)))-x(id)) W/1.4 ]);
		t2 = find( (x>=(S-WL)) & (x<=(S+WR)) )';
		tl = t2(1); tr = t2(length(t2));
		for toremove=t2(:)			% index zone around pos
			wt(toremove) = 0;
			iy(toremove) = y(tl) + ( ( y(tr) - y(tl) ) / ( x(tr) - x(tl) ) * ( x(toremove) - x(tl) ) );
		end	% linear interpolation
%		y=y.*wt;
		plot(x,iy,x,pmax.*y)
		fprintf(1,'\nPeak No %i at X=%f W=(%f,%f) removed. Points %i to %i.', index, S, WL, WR, tl, tr );
	end
end
t2=find(wt ~= 0);
wt=wt(t2);
nx=x(t2);
ny=y(t2);
t1=find( xs < (pos-width) );
t2=find( xs > (pos+width) );

a = xs(t1); b = xs(t2)
c = ys(t1); d = ys(t2);
nx = [ a(:) nx(:) b(:) ];
ny = [ c(:) ny d(:) ];
iy = [ c(:)  iy d(:) ];

fprintf(1,'Length of new data is %i (peaks removed).\n',length(nx));
fprintf(1,'Length of interpolated data is %i (same as initial signal).\n',length(iy));

if (nargout == 1)
	nx=ny;
end

if (size(ys,2)==1)
	nx = nx';
	ny = ny';
	sy = sy';
	iy = iy';
end

