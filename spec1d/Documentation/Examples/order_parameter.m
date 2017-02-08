function order_parameter
%
% This example reads in a set of rocking curves
% used to determine the order parameter of a magnetic transition.
% It fits each curve to a Gaussian, and then fits the resulting
% amplitudes to a power law to determine the order parameter. Notice
% that at no point is it necessary to use a loop structue.
%
% February 2001
% Des McMorrow

% First load data files to extract temperatures PTEM of scans
s1=loads('tasbatch','data/beta000[7:9].dat,X=PTEM,Y=I');
s2=loads('data/beta00[10 12 13 16:35].dat');
stot=[s1 s2];

% Extract temperatures and work out mean
temperature=extract(stot);
temperature=mean(temperature)';

% Now load scan data into two spec1d arrays, and then horizontally concatenate them to
% make one large  spec1d array
s1=loads('data/beta000[7:9].dat,X=OM,Y=I,M=MON');
s2=loads('data/beta00[10 12 13 16:35].dat');
stot=[s1 s2];

% Estimate the peak parameters using the method of moments and plot them
%stat=peakm(stot);
%figure
%plot(temperature,stat(:,4),'ro-')
%xlabel('Temperature')
%ylabel('Integrated Intensity')

% Fit scans to a Gaussian peak, saving results in fitdata

% Either as a loop, so that each fit can be inspected
%for il=1:length(stot)
%   [r,fitdata(il)]=fits(stot(il),'gauss');
%plot(r);pause
%end

% Or as a single operation
[r,fitdata]=fits(stot,'gauss');

% The results are now in the structure fitdata. The parameter values, errors etc are written to an
% ordinary array using
pvals=[fitdata(:).pvals];
evals=[fitdata(:).evals];

% Write fit output to a spec1d object for plotting
sf=spec1d(temperature,pvals(1,:),evals(1,:));
sf=setfield(sf,'x_label','Temperature','y_label','Intensity');
plot(sf)

% Fit order parameter to power law and plot
[sff,fitdata]=fits(sf,'pow',[0.1 78 0.5 0.093],[1 1 1 1]);
plot(sff)
formatpars(fitdata,1);
