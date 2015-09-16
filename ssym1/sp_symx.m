function [newx, centre, peaksav, peakanalysis, rescaletable] = sp_symx(x, y, err, order, option, peakanalysis, centre, cencriteria, peak, appfcn, specfcn);
% [newx, centre, peak, peakanalysis, rescaletable] = 
%    sp_symx(x, y, err, order, option, peakanalysis, centre, cencriteria, peak, appfcn, specfcn);
% Spectra symetrization according to peak positions in signal (x,y) and centre.
% input :
%   order   : polynomial fit level used for final symetrization, can be empty (default is 2).
%   option  : boolean vector for options (default is [ 1 1 0 1 0 ])
%       * symetrization for apparatus (1=yes).
%       * symetrization for signal (1=yes).
%       * global(1)/local(0) fit for apparatus zones.
%       * iterative method (1=yes).
%       * local symetrization (around each peak) (1=yes). First iteration is always global.
%   [centre, peak] : given by sp_gappf, or empty.
%   appfcn  : function name to use to fit apparatus part of spectra (default is gauss).
%   specfcn : function name to use to fit other peaks (default is lorz).

% E.Farhi 04/98 (manuf@ldv.univ-montp2.fr)
% uses : sp_gappf.m, mf_flsqr, mf_simplx, mf_autogs

if ~exist('option') | length(option) < 5
	option = [ 1 1 0 1 0 ];
end
if ~exist('appfcn')
	appfcn = 'gauss';
	option(3) = 0;
end
if ~exist(appfcn)
	appfcn = 'gauss';
	option(3) = 0;
end
if ~exist('specfcn') | ~exist(specfcn)
	specfcn = 'lorz';
end
if ~exist('order')
	order = [];
end

if ~exist('cencriteria')
	cencriteria = 0;
	sp_type = 'user request';
end

fitrout = findobj('tag','mf_FitRoutineFile');
if ~isempty(fitrout)
	fitrout=get(fitrout,'string');
else
	fitrout = 'mf_flsqr';
end

partused = option(1:2);
globapp = option(3);
iterative = option(4);
localsym = option(5);

ly = length(y);
if ~exist('centre') | isempty(centre)
	centre = x(round(ly/2));
end

if ~exist('err') | isempty(err)
	err = sqrt(y-min(y)+min(y(find(y>0))));
end

x=x(:); y=y(:); err = err(:);
newx=x;
xs = x; ys = y;
pmin = []; pmax = []; levely = [];
notexitnow = 1; rescaletable = []; 

limcrit = 2e-2;

% start of iterative part ================================================
disp('Working...');
while (notexitnow)

if ~exist('peak') | ~exist('peakanalysis') | isempty(peak) | isempty(peakanalysis)
	[shutterzones, peak, peakanalysis, levely, sp_type, pmin,pmax] = sp_gappf(x,y,levely,pmin, pmax);
end
if isempty(cencriteria) | (cencriteria > 0) | (cencriteria == 0 & notexitnow > 1)
	[centre, cencriteria] = sp_cent(x,y,peak);
end
if isempty(peakanalysis)
	peakanalysis = ones(size(x));
end

if cencriteria == 0
	limcrit = 5e-2;
end

shutter = 0*peakanalysis;
shutter(find(peakanalysis == 1)) = 1;

peaksav = peak;

% Now fit apparatus or normal peaks with given functions.

% verify symetry
% start from centre
% search nearest peak, opposite nearest peak. 
% if not found pass both 
% else identify type (apparatus) or (signal)
% do autoguess + function autoguess
% fit with required function
% store both positions
%
% fit positions with polynome (deg <=n)

if cencriteria >= 1 & ( ~iterative | (iterative & notexitnow==1) )
	fprintf(1,'Sorry : spectra isn''t symetric enough (%f).\n',cencriteria);
	disp('Normal procees for Symetry Criteria < 1')
	disp('Set Symetry Criteria to zero to force process')
	return;
end

rightpart = find(peak(:,4) >= centre); % indexes in peak matrix of peaks on the right side
leftpart =  find(peak(:,4) < centre);

if isempty(rightpart) | isempty(leftpart) | (sum(partused) == 0)
	disp('Sorry : nothing to do.');
	return
end

% [ index left_width right_width  max_pos Intensity Width ]
i=1;
rescaletable = [];
peak(:,4) = peak(:,4) - centre; % centre is now 0 in peak matrix.
fprintf(1,'Symetrization for X axis on %i peaks.\n',length(rightpart));
totalcrit = 0;
totalfound = [ 0 ];
tp = peak(leftpart,4);
ti = peak(leftpart,5);
tw = peak(leftpart,6);
tc = mean(diff(peak(:,4)))/2; % mean distance on x axis between peaks

while i <= length(rightpart)
	rightpindex = rightpart(i);
	rightindex = peak(rightpindex,1);
	rightpos = peak(rightpindex,4);
	rightwidth = 2*max(peak(rightpindex,2),peak(rightpindex,3)); % in points
	righttype = peakanalysis(rightindex);

	critint = abs(log(abs(ti/peak(rightpindex,5))));
	critwid = abs(log(abs(tw/peak(rightpindex,6))));
	criteria = 0;
	notfound = 1;

	if i == 1 & abs(rightpos) < max(peak(rightpart,4)/100) & abs(centre) > 2*abs(rightpos) % new central peak identified
			rightpos = 0;
	elseif i == 1
			[j,k] = min(abs(x - centre));	% index on x axis
			cp = x(k);		% index on x axis
			rescaletable = [ centre 0 0 k ];	% keep previous central peak
			fprintf(1,'%2i %7.2f {%5.2f} [ central ',0,j,0);
			fprintf(1,'position ] [%7.2f]\n',centre);
%			newx = x-centre;
	end
	fprintf(1,'%2i %7.2f ',i,rightpos);
	
	if rightpos > 0
		leftpindex = find( -tp > rightpos-tc & -tp < rightpos+tc );
		critpos = abs(log(abs(tp/rightpos)));
		if length(leftpindex) >= 1
			[d1,i1] = min(critpos(leftpindex));
			[d2,i2] = min(critint(leftpindex));
			[d3,i3] = min(critwid(leftpindex));
			criteria = d1;
			notfound = 0;
		end
		if length(leftpindex) > 1
			if (d1 > .1) | (d2 > 2) | (d3 > 2)
				criteria = max([d1 d2 d3]);
				notfound = 1;
				leftpindex = leftpindex(1);
			else
				notfound = 0;
				if d1 < 1e-3
					leftpindex = leftpindex(i1);
				elseif (i1 == i2) | (i1 == i3)
					leftpindex = leftpindex(i1);
				elseif (i2 == i3)
					leftpindex = leftpindex(i2);
				else
					leftpindex = [];
				end
				criteria = min([d1 d2 d3]);
			end
		end
		if isempty(leftpindex)
			leftpindex = -1;
			notfound = 1;
			criteria = 0;
		end
	else
		criteria = 0; % central peak
		leftpindex = rightpindex;
		critpos = tp;
		notfound = 0;
	end

	fprintf(1,'{%5.2f} ',criteria);
	if notfound == 0
		totalcrit = totalcrit + criteria ;
	end

	if righttype == 1
		fcn = appfcn;
		fprintf(1,'Apparatus ');
	else
		fcn = specfcn;
		fprintf(1,'Signal    ');
	end
	fprintf(1,'(%s) \t',fcn);

%[tp' ; ti' ; tw' ; critpos' ; critint' ; critwid']
	if notfound | ~isempty(find(totalfound == leftpindex))
		lefttype = -1; % peak not found
		if isempty(find(totalfound == leftpindex))
			fprintf(1,'[no symetry found] ');
		else
			fprintf(1,'[peak re-use     ] ');
		end
	else
		leftindex = peak(leftpindex,1);
		leftpos = peak(leftpindex,4);
		leftwidth = 2*max(peak(leftpindex,2),peak(leftpindex,3));
		lefttype = peakanalysis(leftindex);
		if isempty(totalfound)
			totalfound = leftpindex;
		else
			totalfound = [ totalfound leftpindex ];
		end
	end

	if (lefttype == righttype) & (	(partused(1) & (righttype == 1)) | (leftpindex == rightpindex) | (partused(2) & (righttype == 2)))
		leftfit = []; rightfit = [];
% fit left part --------------------------------------------------
		pin = [];
		selected = 0*x;
		selected(peak(leftpindex,7):peak(leftpindex,8)) = 1;
%fprintf(1,'{%i:%i}',peak(leftpindex,7),peak(leftpindex,8));
		[dummy, name, pnames, pin]=feval(fcn,x,pin,1);	% identify function
		if lefttype == 1
			yf = y.*shutter;	% get all apparatus peaks
		else
			yf = y.*(~shutter).*selected;	% only one peak
		end
		[pin,pmat] = mf_autogs(x, yf, selected, [], peaksav, pnames,[]); % general autoguess in selected region

		[dummy, name, pnames, pin]=feval(fcn,x,pin,yf);	% function autoguess
		if (globapp) & (lefttype == 1)
			selected = 1;
		end

		if ~isempty(fcn) & ~isempty(pmat)
			p = feval(fitrout,x,yf.*selected,1,pin,[],fcn,[]);		% fit
		else
			p = pin;
		end

% now look for parameter 'position' and store.
		n = [];
		pnames = lower(pnames);
		for j=1:size(pnames,1)
			n = [ findstr(pnames(j,:),'energ') findstr(pnames(j,:),'pos') findstr(pnames(j,:),'cent') findstr(pnames(j,:),'freq')  ];
			if ~isempty(n)
%				fprintf(1,'%i ',p(j));
				aleftfit = pin(j);
				if pin(j) & abs(log(abs(p(j)/pin(j)))) < 1e-1 & (p(j) > min(x) & p(j) < max(x))
					leftfit = p(j);
					fprintf(1,'[%7.2f] ',leftfit);
				elseif pin(j)
					leftfit = pin(j);
					fprintf(1,'[%7.2f]*',leftfit);
				else
					leftfit = [];
					fprintf(1,'[?] ',leftfit);
				end
				break;
			end
		end
		if isempty(n)
			fprintf(1,'[not found] ');
		end


% fit right part --------------------------------------------------
		if (rightindex ~= leftindex)
			pin = [];
			selected = 0*x;
%fprintf(1,'{%i:%i}',peak(rightpindex,7),peak(rightpindex,8));
			selected(peak(rightpindex,7):peak(rightpindex,8)) = 1;
			[dummy, name, pnames, pin]=feval(fcn,x,pin,1);	% identify function
			if righttype == 1
				yf = y.*shutter;	% get all apparatus peaks
			else
				yf = y.*(~shutter).*selected;	% only one peak
			end
			[pin,pmat] = mf_autogs(x, yf, selected, [], peaksav, pnames,[]); % general autoguess
			[dummy, name, pnames, pin]=feval(fcn,x,pin,yf);	% function autoguess

			if (globapp) & (righttype == 1)
				selected = 1;
			end

			if ~isempty(fcn) & ~isempty(pmat)
				p = feval(fitrout,x,yf.*selected,1,pin,[],fcn,[]);		% fit
			else
				p = pin;
			end

% now look for parameter 'position' and store.
			pnames = lower(pnames);
			n = [];
			for j=1:size(pnames,1)
				n = [ findstr(pnames(j,:),'energ') findstr(pnames(j,:),'pos') findstr(pnames(j,:),'cent') findstr(pnames(j,:),'freq')  ];
				if ~isempty(n)
					arightfit = pin(j);
%					fprintf(1,'%i ',p(j));
					if pin(j) & abs(log(abs(p(j)/pin(j)))) < 1e-1 & (p(j) > min(x) & p(j) < max(x)) 
						rightfit = p(j);
						fprintf(1,'[%7.2f] ',rightfit);
					elseif pin(j)
						rightfit = pin(j);
						fprintf(1,'[%7.2f]*',rightfit);
					else
						rightfit = [];
						fprintf(1,'[?] ',rightfit);
					end
					break;
				end
			end
			if isempty(n)
				fprintf(1,'[not found] ');
			end
		else
			fprintf(1,'[ centre] ');
			rightfit = leftfit;
%			newx = x-rightfit; 
		end
% now set axisrescale table.
		if ~isempty(rightfit) & ~isempty(leftfit)
			if rightfit ~= leftfit
				rescalepos = (abs(rightfit-centre)+abs(leftfit-centre))/2;
				% we test if it is really symetric...
				if abs(log(abs(rightfit-centre)/abs(leftfit-centre))) < limcrit | ...
					 abs(log(abs(rightfit)/abs(leftfit))) < limcrit
					fprintf(1,'-> [%7.2f] ',rescalepos);
					toadd = [ rightfit rescalepos rightpindex rightindex ; leftfit -rescalepos leftpindex leftindex ];
				else
					rightfit = arightfit; leftfit = aleftfit;
					rescalepos = (abs(rightfit-centre)+abs(leftfit-centre))/2;
					if abs(log(abs(rightfit-centre)/abs(leftfit-centre))) < limcrit | ...
					   abs(log(abs(rightfit)/abs(leftfit))) < limcrit
						fprintf(1,'-> [%7.2f]a',rescalepos);
						toadd = [ rightfit rescalepos rightpindex rightindex ; leftfit -rescalepos leftpindex leftindex ];
					else
						fprintf(1,'-> [~sym:%7.2f]',rescalepos);
						toadd = [];
					end
				end
			else
				centre = rightfit;
				fprintf(1,'-> [   0.00] ');
				[j,k] = min(abs(x - rightfit));	% index on x axis
				toadd = [ rightfit 0 rightpindex k ];
			end
			if isempty(rescaletable)
				rescaletable = toadd;
			else
				rescaletable = [ rescaletable ; toadd ];
			end
		else
			fprintf(1,'- Passed');
		end
	else
		fprintf(1,'- Passed');
	end
	i = i +1;
	fprintf(1,'\n');
end

% now do polynomial fit for mean positions
[dummy, sorti] = sort(rescaletable(:,1));

if isempty(order) 
  order = 2;
end
if order == Inf
	order = size(rescaletable,1)-1; % exact
end
if order >= size(rescaletable,1) | order < 1
	order = size(rescaletable,1); % exact
	order = ceil(order/2);
end

if size(rescaletable,1) <= 3
	if isempty(rescaletable)
		return
	else
		newx = x-rescaletable(1,1);
		disp('Shifting centre');
	end
else
if ~localsym % | notexitnow == 1
	rescaletable = rescaletable(sorti,:);
	fprintf(1,'Global polynomial fit (order %i, %i points) and axis centering\n',order,size(rescaletable,1));
	px = rescaletable(:,1); py = rescaletable(:,2); pe = err(rescaletable(:,4));
	pol = polyfit(px,py,order); % guess
	pol = feval(fitrout,px,py,pe,pol,[],'polynomial',[]);	% fit
	newx = polyval(pol,x);
	peaksav(:,4) = polyval(pol,peaksav(:,4));
else
	if order > 2, order=2; end
	fprintf(1,'Local polynomial fit (order %i, %i points) and axis centering\n',order,size(rescaletable,1));

	% get centre
	xc = rescaletable(1,1); % old centre position

	[j,i] = min(abs(peak(:,4) - xc));	% index in peak matrix
	cp = peak(i,1);		% index on x axis
	cwl = 2*peak(i,2);
	cwr = 2*peak(i,3);
	lastleftindex = cp-cwl;
	lastrightindex = cp+cwr;
	newx = x - xc;	% just translate for centre

	for lefti=3:2:size(rescaletable,1)
		i = rescaletable(lefti,3);
		leftwl = 3*peak(i,2);
		leftwr = 3*peak(i,3);
		leftp  = peak(i,1);
		if lastleftindex < leftp
			disp('A peak is passed ? left')
[ leftp lastleftindex cp ]
		end
		leftindex = min(lastleftindex,(leftp+leftwr));
		leftdom= max(1,(leftp-leftwl)):leftindex;
		xl     = rescaletable(lefti,1);	% old
		kl     = rescaletable(lefti,2);	% new

		righti = lefti-1;
		i = rescaletable(righti,3);
		rightwl = 3*peak(i,2);
		rightwr = 3*peak(i,3);
		rightp  = peak(i,1);
		if lastrightindex > rightp
			disp('A peak is passed ? right')
[ cp lastrightindex rightp ]
		end
		rightindex = max(lastrightindex,(rightp-rightwl));
		rightdom= rightindex:min(length(x),(rightp+rightwr));
		xr      = rescaletable(righti,1);	% old
		kr      = rescaletable(righti,2);	% new

		px = [ xl xc xr ]; py = [ kl 0 kr ];
		pol = polyfit(px,py,order); % get rescaling polynome
		newx(leftdom) = polyval(pol,x(leftdom));
		newx(rightdom) = polyval(pol,x(rightdom));

		if leftindex < lastleftindex
			newx(leftindex:lastleftindex) = linspace(newx(leftindex),newx(lastleftindex),lastleftindex-leftindex+1);
		end
		if rightindex > lastrightindex
			newx(lastrightindex:rightindex) = linspace(newx(lastrightindex),newx(rightindex),rightindex-lastrightindex+1);
		end
		lastleftindex = max(leftdom);
		lastrightindex = min(rightdom);
	end
end
end % if size rescaletable

fprintf(1,'Mean axis shift : %.2f\n',std(x-newx));
fprintf(1,'Symetry quality : %g\n',totalcrit/length(rightpart));

% iteration handling

if ~iterative | (size(rescaletable,1) <= 3)
	notexitnow = 0;
else
	if  (notexitnow < 5) & (std(x-newx)/ly > 1e-4)
		notexitnow = notexitnow + 1;
		fprintf(1,'Iteration %i...\n',notexitnow);
		x = newx;
		peak = []; pmin = []; pmax = []; peakanalysis = [];
	else
		fprintf(1,'Iteration process finished (%i)',notexitnow-1);
		notexitnow = 0;
	end
end

end

% end of iterative part ================================================
