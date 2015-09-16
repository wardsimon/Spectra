function varargout=getfield(s1,varargin)
%
% function varargout=getfield(s1,varargin)
%
% @SPEC1D/GETFIELD gets the fields in a spec1d structure
%
% Example: >>r=gfield(s103,'x_label','x')
%            r is a cell array
%
% DFM 1.3.99
%
if ~isa(s1,'spec1d'); disp('Not a spec1d object'), return, end

%----- Convert s1 to structure and delete non-numeric fields
  
i=1;
s1_cell=struct2cell(s1);  
while i<=length(varargin)

  if ischar(varargin{i})
     if strmatch(varargin{i},fieldnames(s1),'exact'); 
        cell_num=strmatch(varargin{i},fieldnames(s1),'exact');
        varargout(i)={s1_cell{cell_num}};
     end
  end

i=i+1;
end


