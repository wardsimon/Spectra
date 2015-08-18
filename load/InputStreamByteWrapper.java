/*********************************************************************
Java Wrapper to the java.io.InputStream class for Matlab
This class allocates a 4096 bytes buffer, and call read method
The buffer is returned to Matlab

Compile with: javac InputStreamByteWrapper.java

   % Matlab example that reads a file using Java code and writes it
   % back to a temporary file using Matlab code. Finally the contents
   % of the new file is displayed.

   reader = InputStreamByteWrapper; % Default buffer size is 4096 bytes.

   in = java.io.FileInputStream('InputStreamByteWrapper.java');

   bfr = [];
   len = 1;
   while (len > 0)
     len = reader.read(in, 16); % Read 16 bytes at the time (offset=0).
     if (len > 0)
       bfr = [bfr; reader.bfr(1:len)]; % Add bytes to my Matlab buffer.
     end
   end

   close(in);
   clear in, reader;

   disp(bfr');

   tmpfile = tempname;
   fh = fopen(tmpfile, 'wb');
   fwrite(fh, bfr, 'char');
   fclose(fh);

   type(tmpfile);

Origin: comp.soft-sys.matlab discussion started on 2001-01-04 08:10:30 PST
Solution to the problem that Matlab  arguments to Java are "copied 
by value" and not by reference.

Author: Henrik Bengtsson, Lund University, Sweden
*********************************************************************/
import java.io.*;

public class InputStreamByteWrapper {
     public static byte[] bfr = null;

     public InputStreamByteWrapper(int capasity) {
       bfr = new byte[capasity];
     }

     public InputStreamByteWrapper() {
       this(4096);
     }

     public int read(InputStream in, int offset, int length) throws
   IOException {
       return in.read(bfr, offset, length);
     }

     public int read(InputStream in, int length) throws IOException {
       return read(in, 0, length);
     }
   }
