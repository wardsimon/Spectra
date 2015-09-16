function r=setfield(s1,varargin)
%
% function r=setfield(s1,varargin)
%
% @SPEC1D/SFIELD sets the fields in a spec1d structure
%
% Example: >>r=sfield(s103,'x_label','New Label','x',[1 2 3 4 5])
%
% DFM 1.3.99
% %


s = struct(s1);
r = spec1d(s,varargin{:});

% if ~isa(s1,'spec1d'); disp('Not a spec1d object'), return, end
% 
% %----- Covert s1 to structure and delete non-numeric fields
%   
% i=1;
% s1_cell=struct2cell(s1);
% while i<=length(varargin)
% 
%   if ischar(varargin{i})
%      if strmatch(varargin{i},fieldnames(s1),'exact'); 
%         cell_num=strmatch(varargin{i},fieldnames(s1),'exact');
%         if strmatch(class(s1_cell{cell_num}),class(varargin{i+1}),'exact')
%            if  strmatch(class(varargin{i+1}),'double','exact')                   
%                if length(varargin{i+1})==length(s1_cell{1}) | isempty(varargin{i+1})
% %                  disp('Numeric field OK')
%                   s1_cell{cell_num}=varargin{i+1}(:);
%                   i=i+1;
%                else
%                   disp('Error: arrays must be the same size!')
%                end
%            elseif  strmatch(class(varargin{i+1}),'char','exact')                   
% %              disp('Character field OK')
%               s1_cell{cell_num}=varargin{i+1};
%               i=i+1;
%            end
%         else
%            disp('Error: Classes must be the same!')
%         end
%      else
% %        disp('Error: field does not exist')
%      end
%    end
% i=i+1;
% end
% 
% %--- Create return spectrum
% 
% r.x=s1_cell{1};
% r.y=s1_cell{2};
% r.e=s1_cell{3};
% r.x_label=s1_cell{4};
% r.y_label=s1_cell{5};
% r.datafile=s1_cell{6};
% r.yfit=s1_cell{7};
% 
% r=spec1d(r);