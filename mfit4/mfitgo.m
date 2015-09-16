function mfitgo(x, y, err, xlab, ylab, name)
%
% MFIT function mfitgo
%     Interface to mfit for other programs - e.g. mview
%     MZ 24.6.96
%

%----- Start mfit if it isn't already running ------------------------
if (isempty(findobj('tag','mf_ControlWindow')))
  mfit;
end

%----- Check that the number of parameters are right --------------------------
if nargin<3
        errordlg('Mfitgo called with bad number of parameters');
        return
end


%----- Check that vector sizes match and are non-zero -------------------------
x=x(:); y=y(:); err=err(:);
if ((length(x)~=length(y)) | (length(x)~=length(err)) | length(x)==0)
  errordlg('Mfitgo called with bad parameters');
end

%----- Create axis labels -----------------------------------------------------
if nargin <4, xlab=''; end
if isempty(xlab)
  xlab='x';
end
if nargin <5, ylab=''; end
if isempty(ylab)
  ylab='y';
end
if nargin <6, name=''; end
if isempty(name)
  name = '';
end

%------------ Sort data ---------------------------
[x, i]=sort(x);
y=y(i);
err=err(i);

%------------ Eliminate data with zero error -------
i=find(err==0);
if ~isempty(i)
  mf_msg('Data with zero error adjusted');
  if ~isempty(find(err))
    err(i)=min(err(find(err)))/10;
  else
    err = ones(size(y))*mean(y)/100;
  end
end

%------------ Make data window ---------------------

disp(['Mfit Loading ' name ' : ' ylab ' vs ' xlab ]);
mf_msg('Plotting data');
hmf_data=mf_dwin(xlab, ylab);figure(hmf_data);

%---------- Attach data to userdata ------------------
set(hmf_data,'userdata',[x y err ones(size(y))]);

%----------------- Set limits  --------------------
yrange=max(y)-min(y);
if yrange==0
   yrange=yrange+max([1e-3*mean(y) 1e-6]);
end
xrange=max(x)-min(x);
if xrange==0
   xrange=xrange+max([1e-3*mean(x) 1e-6]);
end

set(gca,'Xlim',[min(x)-0.01*xrange max(x)+0.01*xrange]);
set(gca,'Ylim',[min(y)-0.1*yrange max(y)+0.1*yrange]);

%----------- Do the plot ------------------------------
%mf_uplot;
if ~isempty(name)
  title(name);
  set(findobj('Tag','mf_DataFile'),'string', name);
end
mf_gdata('noload');
mf_dwin(xlab, ylab);



