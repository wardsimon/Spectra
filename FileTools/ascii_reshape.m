function ascii_reshape(filename,string1,string2)
% ascii_reshape : Change characters in a file
%Syntax: ascii_reshape(filename,string1,string2)
% This function changes characters of file 'filename' matching 'string1'
% into 'string2' characters.
% WARNING : The file will be modified.
% Example: ('name',',;','. ') wil change ',' into '.' and ';' into spaces.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Change characters in a file

% uses : none
% E.Farhi 03/96

if (nargin ~= 3)
        error('usage : ascii_reshape(filename,string1,string2)');
end

if (length(string1) ~= length(string2))
        error('strings must be of same length.');
end

[lsal, fnum] = system([ 'ls -l ' filename ]);

if (fnum ~= 0)
        error('data loading : file not found.');
end

% Loading file

fnum = fopen(filename,'r');
[str, fl] = fread(fnum, Inf, 'char');   % all file
fclose(fnum);

implicit_str_to_num_ok = 1;

str = str + 0;          % casts char to int.
string1 = string1 + 0;
string2 = string2 + 0;

for i=1:length(string1)
  pos = find( (str - string1(i)) == 0);
  for j=pos'
      str(j) = string2(i);
  end
  fprintf(1,'\n%i characters "%c" changed into "%c"', length(pos), setstr(string1(i)),  setstr(string2(i)));
end

% str = setstr(str);

fnum = fopen(filename,'w');
fl = fwrite(fnum, str, 'char');
fclose(fnum);
fprintf(1,'\n%i bytes written. \n', fl);


%end
