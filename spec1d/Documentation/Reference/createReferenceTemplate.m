function createReferenceTemplate(fun)

out = sprintf('reference_%s.m',fun);
fid = fopen(out,'w');

if fid <=0
    error('Could not write file')
end

fprintf(fid,'%%%% @spec1d/%s\n',fun);
fprintf(fid,'%% This is the reference documentation for the function @spec1d/%s\n',fun);
fprintf(fid,'%%\n');
fprintf(fid,'%% This function \n');
fprintf(fid,'%%%%\n');
fprintf(fid,'\n');
fprintf(fid,'%%%% Syntax\n');
fprintf(fid,'%%\n'); 
fprintf(fid,'%%    s_out = %s(s)\n',fun);
fprintf(fid,'\n');
fprintf(fid,'%%%% Inputs\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% * _s_ - Single or vector of spec1d objects\n'); 
fprintf(fid,'%%\n'); 
fprintf(fid,'\n');
fprintf(fid,'%%%% Outputs\n');
fprintf(fid,'%%\n'); 
fprintf(fid,'%% * _s_out_ - Spec1d object of the same size as the input array\n');
fprintf(fid,'\n');
fprintf(fid,'%%%% Example\n');
fprintf(fid,'%% This is an example on using @spec1d/%s\n',fun);
fprintf(fid,'%%\n');
fprintf(fid,'%% <html><h3>Example 1</h3></html>\n');
fprintf(fid,'%%\n');
fprintf(fid,'\n');
fprintf(fid,'s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);\n');
fprintf(fid,'s_1 = %s(s);\n',fun);
fprintf(fid,'\n');
fprintf(fid,'figure\n');
fprintf(fid,'plot(s,s_1)\n');
fprintf(fid,'legend({''s'',''%s(s)''})\n',fun);
fprintf(fid,'\n');
fprintf(fid,'%%%% See Also\n');
fprintf(fid,'%% <html><a href="{{ site.url }}/@spec1d.plus/index.html">Plus</a>, <a href="{{ site.url }}/@spec1d.sum/index.html">Sum</a></html>');

fclose(fid);

edit(out)