function tf = ishg2(fig)
try
    tf = ~graphicsversion(fig, 'handlegraphics');
catch
    tf = false;
end