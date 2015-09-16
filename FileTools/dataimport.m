function [data,numtable,choosed] = dataimport(file,matlabgui)
% dataimport : General data importation
%Syntax: [data,numtable,choosed] = dataimport(file)
% Forms a data matrix containing all the fields you choosed in a text data file

% Author:  EF <manuf@ldv.univ-montp2.fr> 27.06.97
% Description:  General load routine. Get all data.

% needs looktxt program

if exist('findobj')
	matlabgui = 1;
else
	matlabgui = 0;
end

%--------- Open data file ---------------------------------

fid=fopen(file,'r');
if (fid<0)
   error('File not found');
   return;
end

filestr = fread (fid);
filestr = setstr(filestr');
fclose(fid);

%------  looktxt analysis ---------------------------------

cmd = [ 'looktxt -F="lktmp001" -t -r -a -f -v ' file ];
system( cmd );

if matlabgui
	patsav = path;
	addpath(pwd);
	path(path);
end
eval('lktmp001');
if matlabgui
	path(patsav);
	path(path);
	eval('delete lktmp001.m');
else
	system('rm lktmp001.m');
end

% now table is 'table'
% [ id type start_pos end_pos rows columns ]
% with type = 1 (num) | 2 (text) | 64 (comment)
numtable = table(find(table(:,2) == 1),:);
[n,c] = sort(numtable(:,5).*numtable(:,6));
numtable = numtable(c,:);
[maxel,maxid] = max(numtable(:,5).*numtable(:,6));
[n,c] = size(numtable);
if (n > 30)
        fprintf(1,'Warn : Too many (%i) fields... \n',n);
        fprintf(1,'Press Ctrl-C to abort if needed\n\n');
end

endheaderpos = min([ numtable(maxid,3) numtable(maxid,4) ]);
if (endheaderpos > 500)
        endheaderpos = min([ numtable(1,3) numtable(1,4) ]);
end

header = filestr(max(1,endheaderpos-500):endheaderpos);
fprintf(1,'Header : ...\n%s\n',header);

choosed = [];
if (n>1)
  fprintf(1,'Found %i fields\n',n);

        for i=1:n
                fprintf(1,'%i - n%i = (%ix%i)\n',i,numtable(i,1),numtable(i,5),numtable(i,6));
        end
        choosed = input ('Enter field numbers you wan"t to import : ','s');
        choosed = str2num(choosed);

else
        fprintf(1,'Found one field... ok\n');
        choosed = 1;
end

if (isempty(choosed))
                disp('Nothing imported.');
        return;
end

data = [];
n = length(choosed);
fprintf(1,'Getting %i fields : ',n);

nr = sum(numtable(choosed,5));
nc = max(numtable(choosed,6));
data = zeros(nr,nc);
r=1;
for i=1:n
        nr = numtable(choosed(i),5);
        nc = numtable(choosed(i),6);
        ni = ['n' num2str(numtable(choosed(i),1))];
        fprintf(1,'%s  ',ni);
        eval([ 'nf = ' ni ';' ]);
        [nr,nc] = size(nf);
        data((r:(r+nr-1)),1:nc) = nf;
        r=r+nr;
end
[nr,ncol] = size(data);
fprintf(1,'\nData matrix has %ix%i elements.\n',nr,ncol);
n = min(10,ncol);
for i=1:min(2,nr)
        fprintf(1,'%g ',data(i,1:n));
        fprintf(1,'\n');
end
fprintf(1,'    ...\n');
for i=max(1,nr-1):nr
        fprintf(1,'%g ',data(i,1:n));
        fprintf(1,'\n');
end


