function [x, y, err, xlab, ylab,monitor] = frombase(file,y,err,name2)
% MFIT/MView direct entry : Get new data from user
% just fill in the x,y and err fields
% with your expressions
%
% You can use your own Matlab variables
% Use 'mf_uplot(p,dp)' to send parameters from Matlab to Mfit
% EF 03.07.97
% [x, y, err, xlab, ylab,monitor] = frombase(file)
% or frombase('mfit') or frombase('base') or frombase(x,y,err)

% vars : x y err i xrange yrange

% direct equivalent using global vars for :
%[x y err xlab ylab]=feval(loadfun, [datadir datafile]);

% first save global variables with same name as the ones now used...

if nargin < 4, name2=''; end
if nargin < 3, err=''; end
if nargin < 2, y=''; end
if nargin < 1, file=''; end
x=[];

if ~isstr(file)
  x = file;
  if ~isempty(name2)
    file = name2;
  else
    file = 'LOCAL :';
  end
end

xlab='';ylab=''; p=[]; fit=[];

if isempty(file)
  file = 'BASE :';
end

if strcmp(file,'mfit')
  name = 'MFIT : Modify data';
else
  name = [ file ' Enter data' ];
end

% get and save base data base -----------------------------------

clx = 0; cly = 0; clerr = 0; clp = 0; clfit = 0;
tmp = 'save tmptomf.mat ';

if (isempty(x) & evalin('base','exist(''x'')') == 1)
  tmp = [ tmp 'x ' ];
  x=evalin('base','x');
  clx = 1;
end

if (isempty(y) & evalin('base','exist(''y'')') == 1)
  tmp = [ tmp 'y ' ];
  y=evalin('base','y');
  cly = 1;
end

if (isempty(err) & evalin('base','exist(''err'')') == 1)
  tmp = [ tmp 'err ' ];
  err=evalin('base','err');
  clerr = 1;
end
if (evalin('base','exist(''p'')') == 1)
  tmp = [ tmp 'p ' ];
  p=evalin('base','p');
  clp = 1;
end

if (evalin('base','exist(''fit'')') == 1)
  tmp = [ tmp 'fit ' ];
  fit=evalin('base','fit');
  clfit = 1;
end

if (clx+cly+clerr+clp+clfit > 0)
  evalin('base',tmp);
end

% get mfit data if needed -----------------------------------

if strcmp(file,'mfit')
  [xf,yf,errf,selected,fitf,pf,dp,fixed]=fromfit('silent');
  if ~isempty(xf) | ~isempty(pf)
    fprintf(1,'MFit data collected : ');
  end
  if ~isempty(xf), fprintf(1,'x '); x=xf; end
  if ~isempty(yf),  fprintf(1,'y '); y=yf; end
  if ~isempty(errf), fprintf(1,'err '); err=errf; end
  if ~isempty(pf), fprintf(1,'p '); p=pf; end
  if ~isempty(fitf), fprintf(1,'fit '); fit=fitf; end
  fprintf(1,'\n');
end
disp('Matlab x,y,err user variables can be used in expressions')

% now do evals ------------------------

assignin('base','x', x);
assignin('base','y', y);
assignin('base','err', err);
assignin('base','p', p);
assignin('base','fit', fit);

[ xrange, yrange, i]=mf_exprdg(name); % ask user

if ~isempty(xrange)
  x=evalin('base', xrange, '[]');
end

if ~isempty(yrange)
  y=evalin('base', yrange, '[]');
end
y=y(:);
assignin('base','y', y);
if isempty(x), x=1:length(y); end
x=x(:);
assignin('base','x', x);

if ~isempty(i)
   err=evalin('base', i,'[]');
end

if (isempty( err) | all( err == 0))
       err=1e-3*max(abs( y));       % no error bars - so equal weights
end

if (length(err) == 1)
  err = err * ones(size(y));
end

err= err(:);

%------------ Exit if load unsuccessful -------------

evalin('base','clear x y err p fit xlab ylab');

if (clx+cly+clerr+clp+clfit > 0)
  evalin('base','load tmptomf.mat');
  delete('tmptomf.mat');
end

i=[ length(x); length(y); length(err) ];
if any(find(i == 0))
  disp('nothing loaded');
  return;
end

%------------ Sort data ---------------------------
[ x,  i]=sort( x);
 y= y( i);

%------------ Eliminate data with zero error -------
 i=find( err==0);
if ~isempty( i)
  disp('Data with zero error eliminated');
   x( i)=[];
   y( i)=[];
   err( i)=[];
  if isempty( x)
     error('Data has zero error');
   end
end


%----------- restore base variables ------------------------------


if strcmp(file,'base')
  disp(sprintf('*Data Loaded. %d points.\n', length( x)));
end
if strcmp(file,'mfit')
  disp('Use ''mf_upars(p,dp)'' to send parameters from Matlab to Mfit');
end
xlab=xrange; ylab = yrange;
monitor = ones(size(y));
