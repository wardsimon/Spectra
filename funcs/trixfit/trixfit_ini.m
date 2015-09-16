function varargout = trixfit_ini(initialise)
%
% TRIXFIT function to initialise all parameters
% 
% xsec:               cross-section file
% method:             rc_cnmat for Cooper-Nathans, rc_popma for Popovici
% parameter_file:     rescal parameter file
% configuration_file: configuration file for Popovici, can be specified as [] if 
%                     rc_cnmat is selected
%

global global_trixfit
global first_call

error_status=[];

global_trixfit.monitor_flag=initialise.monitor_flag;                   
global_trixfit.monte_carlo_samples=initialise.monte_carlo_samples;            

global_trixfit.bkgd_file=initialise.bkgd_file;
global_trixfit.pnam_file=initialise.pnam_file;
global_trixfit.corr_file=initialise.corr_file;
if isfield(initialise,'parallel')
    global_trixfit.parallel = initialise.parallel;
else
    global_trixfit.parallel = 0;
end
if isfield(initialise,'low_pass_filter')
    global_trixfit.low_pass_filter = initialise.low_pass_filter;
else
    global_trixfit.low_pass_filter = 1;
end
if isfield(initialise,'oversample')
    global_trixfit.oversample = initialise.oversample;
else
    global_trixfit.oversample = 0;
end

global_trixfit.R_Do = 1;

parameter_file=initialise.rescal_pars;
configuration_file=initialise.popovici_pars;
xsec=initialise.xsec_file;
method=initialise.resolution_method;

%----- Set defaults for input

if nargin==0
   disp('Must specify cross-section')   
   error_status=1;
   return
end

%----- Set initial value of flag used in triplefit

first_call=1;

%----- Initialise global parameters

pcn=[]; ppop=[];

%----- Read in CN parameters
fid=fopen(parameter_file);
if fid<0
   warning('Cooper-Nathans parameter file not in path')
   error_status=1;
   return
end
pcn=parameter_read(fid);
if length(pcn)~=42
   warning('Cooper-Nathans requires 42 parameters!')
   error_status=1;
   return
end

%----- Read in Popovici parameters

if isempty(configuration_file) & strcmp(method,'rc_popma')
   warning('Configuration file must be specified for Popovici method')
   error_status=1;
   return
end

if ~isempty(configuration_file) & ~strcmp(method,'rc_cnmat')
    fid=fopen(configuration_file);
    if fid<0
        warning('Popovici parameter file not in path')
        error_status=1;
        return
    end
    ppop=parameter_read(fid);
    new_rescal=strfind(which('rescal.m'),'rescal6');
    
    if length(ppop) == 30
        if ~new_rescal
            warning('The input file has a description of the monitor. This is now neglected as you are running rescal5!')
            ppop(28:30)=[];
        end
    elseif length(ppop) == 27
        if new_rescal
            warning('The input needs a description of the monitor. This is as you are running rescal6!')
        end
    else
        warning('Popovici method requires 27 or 30 parameters. Not %i',length(ppop)) % including the monitor parameters, the .cfg file has 30 parameters
        error_status=1;
        return
    end
    
else
    ppop=zeros(27,1);
end

%----- Define global variables

%----- Check to see if cross-section exists

check_method=exist(xsec);
if sum(check_method==[2 3])<1
   warning('Cross-section not in path')
   error_status=1;
   return
end

%---- Check to see if method exists

check_method=exist(method);
if check_method~=2
   warning('Resolution calculation method not in path')
   error_status=1;
   return
end

global_trixfit.xsec_file=xsec;
global_trixfit.resolution_method=method;
global_trixfit.pres=[pcn; ppop];
global_trixfit.first_call=1;

if ~isempty(error_status)
    warning('Trixfit initialisation error.....')
end

if nargout == 1 
    varargout{1} = error_status;
end

%=======================================================
function [p]=parameter_read(fid)
%
%
%

%-------------- Initialize arrays-----------------

data=[];		
header='';
text=fgetl(fid);

%----------------- Load data-----------------------

while (text>0)
   [temp count]=sscanf(text,'%f');			
   if isempty(temp) 
      header=[header text];
   else
      if (count==size(data,2) | isempty(data))
         data=[data; temp'];
      end
   end
      text=fgetl(fid);
end
fclose(fid);

p=data(:,1);
