/* --------------------- Program looktxt.c ---------------------- */

/* 
Usage : looktxt [options] filename

 Action: Search and export numerics in a text/ascii file.
 The programs looks into your file some numeric fields (and optionally
 characters). Each identified numeric field is then named and exported
 into an output filename. The output file can be a .m (Matlab), .oct
 (Octave) or .txt file, with automatic or specified user name. For a
 special type of data, you can specify point, separator, end of line
 and comment conventions. A 'table' field can be exported, containing
 informations about all other found fields (use looktxt -h for help). 
 (c) I.L.L. Grenoble, France.
 URL http://www.ill.fr
*/  

/* content: C language
* compile with : gcc -O looktxt.c -o looktxt */

/* History:
* 0.86  (04/11/97) not effective
* 0.87  (26/03/99) works quite fine. Some bugs.
* 0.88  (09/04/99) improvements and 'table' output
* 0.89  (02/07/99) corrected grouping error for isolated numerics 
* 0.89a (27/03/00) multi plateform handling
* 0.90  (03/07/00) debug mode ok, no more lktmp00 file
* 0.91  (26/07/00) new options -S (struct) -H (num header)
* 0.93  (21/08/01) -T, filename in file
*/

#define author  "Farhi E. 01/97 <farhi@ill.fr>"
#define date    "22/11/01"
#ifndef version
#define version "0.9.4"
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef WIN32
#define LK_PATHSEP_C '\\'
#define LK_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define LK_PATHSEP_C ':'
#define LK_PATHSEP_S ":"
#else  /* !MAC */
#define LK_PATHSEP_C '/'
#define LK_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */

/* Usage :
* ARGS : see help in main
* OUT : matrix of fields sent to stdout (printed to terminal)
can be used after to extract particular data.

 * Creates an '.m', '.oct' or simple numeric output file
 
  * Example :
  looktxt foo.txt -p="." -s='\t\v,;' -c='#'
  
   'looktxt -h' to get help.
*/

/* remove comments if your compiler does not auto include some libraries */
#include <stdio.h>
/* #include <stdlib.h> */
#include <string.h>
/* #include <ctype.h> */
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

/* -- Declaration section ----------------------------------------- */


#define Bnumber      1
#define Balpha       2
#define Bpoint       4
#define Beol         8
#define Bexp        16
#define Bsign       32
#define Bcomment    64
#define Bseparator 128

#define OUT_OCT      1
#define OUT_M        2
#define OUT_TXT      3

/* number format : [+-][0-9][.][0-9][eE][+-][0-9] */

/* event types                           {0,a,.,l,e,-,c,s} */
const long     all         =  0xFF;   /* {1,1,1,1,1,1,1,1}; */
const long     none        =  0x00;   /* {0,0,0,0,0,0,0,0}; */

/* events for num search                 {n,a,p,l,e,s,c,s} */
const long needforstartnum =  Bnumber + Bpoint + Bsign;
const long needafternumber =  Bnumber + Bpoint + Bexp + Beol + Bseparator;
const long needaftersign   =  Bnumber + Bpoint;
const long needafterexp    =  Bnumber + Bsign;
const long needafterpoint  =  Bnumber + Bexp   + Beol + Bseparator; /* but Bnumber only, if num started with point ! */
const long needaftereol    =  Bnumber + Bpoint + Bsign + Beol + Bseparator;
const long needaftersep    =  Bnumber + Bpoint + Bsign + Beol + Bseparator; 


signed char verbose  = 0;
char passascii= 1;
char groupnum = 0;
char noroot   = 0;
char outtype  = OUT_M;
char forcefile= 0;
long line     = 0;
char tablextr = 0;
char username = 0;
char numheader= 0;
char usestruct= 0;
char *userfile;
char *filename;

#define number    "0123456789"
#define alpha     " !\"#$&'()*+,-./123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
#define Cpoint     "."
#define Ceol       "\n\f"
#define exp       "eE"
#define sign      "+-"
#define Ccomment   "#%"
#define Cseparator "\t\v\r,; ()[]{}=|:<>&\"/"

char *point, *comment, *separator, *eol;


char *make_esc_char(char *s)
{  /* make escape char such as \[abfnrtv\].... */
  int i=0;  /* in s */
  int j=0;  /* in output */
  char c1 = '\0';
  char c2 = '\0';
  char *sout;
  
  sout = (char *)malloc(strlen(s)+1);
  do
  {
    c1 = s[i];
    c2 = s[++i];
    if (c1 == '\\')
    {
      switch (c2)
      {
      case 'a' : c1 = '\a'; i++; break;
      case 'b' : c1 = '\b'; i++; break;
      case 'f' : c1 = '\f'; i++; break;
      case 'v' : c1 = '\v'; i++; break;
      case 't' : c1 = '\t'; i++; break;
      case 'n' : c1 = '\n'; i++; break;
      case 'r' : c1 = '\r'; i++; break;
      case '\\': c1 = '\\'; i++; break;
      }
    }
    sout[j++] = c1;
    if (i >255)
    {
      printf("looktxt: warn : make esc char : string too long\n");
      free(sout);
      return(NULL);
    }
  }
  while ((c2 != '\0') && (c1 != '\0'));
  sout[j] = '\0';
  return(sout);
}

char *undo_esc_char(char *s)
{  /* undo escape char such as \[abfnrtv\].... */
  int i=0;  /* in s */
  long j=0;
  char c2 = '\0';
  char *sout;
  
  sout = (char *)malloc(strlen(s)*5);
  sout[0]='\0';
  do
  {
    c2 = s[i++];
    switch (c2)
    {
    case '\a' : sout = (char *)strcat(sout,"\\a"); break;
    case '\b' : sout = (char *)strcat(sout,"\\b"); break;
    case '\f' : sout = (char *)strcat(sout,"\\f"); break;
    case '\v' : sout = (char *)strcat(sout,"\\v"); break;
    case '\t' : sout = (char *)strcat(sout,"\\t"); break;
    case '\n' : sout = (char *)strcat(sout,"\\n"); break;
    case '\r' : sout = (char *)strcat(sout,"\\r"); break;
    default   : j=strlen(sout); sout[j] = c2; sout[j+1]='\0'; break;
    }
    if (i > 255)
    {
      printf("looktxt: warn : undo esc char : string too long\n");
      free(sout);
      return(NULL);
    }
  }
  while (c2 != '\0');
  return(sout);
}


int ischr(char c, char *category)
{
  return (strchr( category   ,c) != NULL);
}

void printfield(char *filestr, long pos, long filesize, long size)
{
  int i;
  
  for (i = 0; ((pos+i < filesize) && (i < 20) && (i <= size)); i++)
    if (isprint(filestr[pos+i])) printf("%c",filestr[pos+i]);
    else printf(" ");
}

long scanfile(char *filename)
{
  FILE *fnum;
  FILE *fout;
  char c;
  char *fileout;
  char *filestr;
  char *filestrcp;
  char tablestr[1024];
  long *filetable;
  char *filenamecp;
  struct stat stfile;
  long filesize;
  char buffer[4096];
  
  int  last_is     = 0; /* type of preceeding char */
  long last_eolpos = 0; /* last EOL position */
  long last_seppos = 0; /* Separator */
  int  is          = 0; /* current char type */
  int  found       = 0; /* is that char expected ? */
  long pos         = 0; /* current pos */
  long rows        = 0; /* number of rows in matrix */
  long columns     = 0; /* number of columns in matrix */
  long last_columns= 0; /* to test if number of columns has changed since last line */
  long startcmtpos = 0; /* comment field starting pos */
  long startnumpos = 0; /* num field start pos */
  long startcharpos= 0; /* char field start pos */
  long endcmtpos   = 0; /* comment field end pos */
  long endnumpos   = 0; /* num field end pos */
  long endcharpos  = 0; /* char field end pos */
  int  possiblenum = 0; /* we migth be in a number ... */
  int  possiblecmt = 0; /* we are in a comment */
  int  fieldend    = 0; /* detect end of field */
  long id          = 0; /* field number */
  long i           = 0; /* index */
  int  need        = 0; /* what is to be expected for next char */
  char inpoint     = 0; /* if we are after a point */
  char inexp       = 0; /* if we are after an exp */
  long prod        = 0; /* cols * rows */
  long numextr     = 0;
  
  long maxfield    = 0; /* dimensions of biggest matrix */
  long maxprod     = 0;
  
  long filetablesize = 1024;
  long filetablecontent = 0;
  time_t ltime;
  
  
  
  /* file section ----------------------------------------------------- */  
  
  /* open, read then close input file */  
  fnum = fopen(filename,"r");
  if (fnum == NULL)
  {
    printf("looktxt: Fatal : unable to open file %s\n",filename);
    return(-1);
  }
  
  stat(filename,&stfile);
  filesize = stfile.st_size;
  filestr = (char *)malloc(filesize+64);
  if (filestr == NULL)
  {
    printf("looktxt: Fatal : unable to malloc %i bytes for filestr (%s)\n",filesize,filename);
    return(-1);
  }
  if (!fread(filestr, 1, filesize, fnum))
  {
    printf("looktxt: Fatal : unable to read file %s\n",filename);
    return(-1);
  }
  fclose(fnum);
  filestr[filesize] = '\n';
  filestr[filesize+1] = '\0';
  
  /* make a copy of input string, and manage file names */
  
  filestrcp = (char *)malloc(filesize+64);
  if (filestrcp == NULL)
  {
    printf("looktxt: Fatal : unable to malloc %i bytes for filestrcp (%s)\n",filesize,filename);
    return(-1);
  }
  strcpy(filestrcp, filestr);
  filenamecp = (char *)malloc(strlen(filename)+64);
  strcpy(filenamecp, filename);
  for (i=0; i < (long)strlen(filename); i++)
  { if ( (filename[i] > 122) || (filename[i] < 32) || (ischr(filename[i], "!\"#$%&'()*+,-.:;<=>?@[\\]^_`/") && (filename[i] !=  LK_PATHSEP_C)) ) filename[i] = '_'; }
  fileout = (char *)malloc(strlen(filename)+64);
  if (fileout == NULL)
  {
    printf("looktxt: Fatal : unable to malloc %i bytes for fileout %s\n",strlen(filename),filename);
    return(-1);
  }
  
  if (username)
    strcpy(fileout, userfile);
  else
  {
    strcpy(fileout, filename);
    strcpy(userfile, filename);
  }
  if ((fileout[0] >= '0') && (fileout[0] <= '9'))
    sprintf(fileout, "L%s", userfile);
  if ((filename[0] >= '0') && (filename[0] <= '9'))
  {
    strcpy(userfile, filename);
    sprintf(filename, "L%s", userfile);
  }
  if (outtype == OUT_OCT) fileout = (char *)strcat(fileout,".oct");
  if (outtype == OUT_M) fileout = (char *)strcat(fileout,".m");
  if (outtype == OUT_TXT) fileout = (char *)strcat(fileout,".txt");
  
  /* open output file and table file */
  
  if (!stat(fileout,&stfile) && !forcefile)
  {
    printf("looktxt: Info : %s already exists. Use '-f' to force file writing.\n",fileout);
    return(-1);
  }
  if (verbose >= 1) printf("Creating output file %s\n",fileout);
  fout = fopen(fileout,"w+");
  if (fout == NULL)
  {
    printf("looktxt: Fatal : unable to create file %s\n",fileout);
    return(-1);
  }
  filetable = (long *)malloc(6*sizeof(long)*filetablesize);
  if (filetable == NULL)
  {
    printf("looktxt: Warn : unable to malloc %i bytes for filetable. No table.\n",6*sizeof(long)*filetablesize);
    tablextr = 0;
  }
  
  /* get root name for field names */
  
  if (strrchr(filename, LK_PATHSEP_C) != NULL) 
  { 
    filename = strrchr(filename,LK_PATHSEP_C)+1; 
  }
  
  /* comments in out file */
  buffer[0] = '\0';
  if ((outtype == OUT_OCT) || (outtype == OUT_M))
  {
    time(&ltime );
    sprintf(buffer,               "%% This file was processed by looktxt v" version "\n");
    sprintf(buffer+strlen(buffer),"%% File:    %s (from %s)\n",fileout,filenamecp);
    sprintf(buffer+strlen(buffer),"%% Date:    %s\n", ctime( &ltime ));
    strcat(buffer               ,"% Command: looktxt ");
    if (!passascii)    strcat(buffer," -a");
    if (groupnum)      strcat(buffer," -g");
    if (noroot)        strcat(buffer," -r");
    if (outtype == OUT_OCT) strcat(buffer," -o=oct");
    if (outtype == OUT_M)   strcat(buffer," -o=m");
    if (outtype == OUT_TXT) strcat(buffer," -o=txt");
    if (forcefile)     strcat(buffer," -f");
    if (tablextr == 1) strcat(buffer," -t");
    if (tablextr == 2) strcat(buffer," -T");
    if (usestruct)     strcat(buffer," -S");
    if (numheader)     strcat(buffer," -H");
    if (username) sprintf(buffer+strlen(buffer)," -F=\"%s\"",userfile);
    sprintf(buffer+strlen(buffer)," %s\n", filenamecp);
    
    strcat(buffer,"% Fields 'nXX' indicate numbers, and 'sXX' are for strings.\n");
  }
  if (!noroot && (strlen(filename) > 25)) 
  {  
    sprintf(filename, "var%i", (long)time(NULL));
    printf("looktxt : Warn : filename too long, using private ID for fields : %s\n",filename);
    if ((outtype == OUT_OCT) || (outtype == OUT_M))
      sprintf(buffer+strlen(buffer),"%% Private ID for fields is %s\n",filename);
  }
  if (outtype == OUT_OCT)
    sprintf(buffer+strlen(buffer),"%% Type in GNU Octave 'octave> load %s.oct' to import data\n",fileout);
  if (outtype == OUT_M)
    sprintf(buffer+strlen(buffer),"%% Type in Matlab 'matlab> %s;' or eval('%s'); to import data\n",fileout, fileout);
  if ((tablextr) && ((outtype == OUT_OCT) || (outtype == OUT_M)))
  {
    strcat(buffer,"%% A 'table' field description is in ");
    if (!noroot) sprintf(buffer+strlen(buffer),"%s_",filename);
    if (usestruct && ((outtype == OUT_M))) strcat(buffer,"s.");
    strcat(buffer,"table with content : \n");
    strcat(buffer,"%% table = [ id type start_pos end_pos rows columns ]\n");
    strcat(buffer,"%% with type = 1 (num) | 2 (text) | 64 (comment)\n");
  }
  if ((outtype == OUT_OCT) || (outtype == OUT_M)) fprintf(fout,"%s", buffer);
  
  /* filename field in out file */
  strcpy(buffer, "");
  if (outtype == OUT_M)
  {
    if (!noroot)   sprintf(buffer,"%s_",filename);
    if (usestruct) sprintf(buffer+strlen(buffer),"s.");
    sprintf(buffer+strlen(buffer),"filename = '");
    sprintf(buffer+strlen(buffer),"%s", filenamecp);
    fprintf(fout,"%s';\n", buffer);
  }
  if (outtype == OUT_OCT)
  {  
    strcpy(buffer,"# name: ");
    if (!noroot) sprintf(buffer+strlen(buffer),"%s_",filename);
    sprintf(buffer+strlen(buffer),"filename\n");
    sprintf(buffer+strlen(buffer),"# type: string\n");
    sprintf(buffer+strlen(buffer),"# length: %i\n", strlen(filenamecp));
    sprintf(buffer+strlen(buffer),"%s", filenamecp);
    fprintf(fout,"%s\n", buffer);
  }
  
  if (verbose >= 1) 
  {
    printf("Scanning...\n");
    printf("   id type     start        end       rows    columns      line\n");
  }
  
  /* init scan ---- */
  
  rows         = 0;
  columns      = 0;
  last_columns = 0;
  is           = Beol;
  startcmtpos  = filesize;
  endcmtpos    = -1;
  startcharpos = 0; endcharpos = 0;
  fieldend     = 0;
  
  
  /* scan --------- */
  
  do
  {
    last_is = is;
    
    c=filestr[pos];
    
    if (pos > filesize) /* end of file reached : real exit */
    {
      last_eolpos = filesize;
      last_seppos = filesize;
      last_is     = Beol;
      possiblecmt = 0;
      need = 0; /* generates end of field : end of line */
      c = '\n';
    }
    if (pos == filesize) /* end of file reached */
    {
      c = '\n';
    }  
    
    is =   Bnumber    * ischr(c, number   )
      + Balpha     * ischr(c, alpha    )
      + Bpoint     * ischr(c, point    )
      + Beol       * ischr(c, eol      )
      + Bexp       *(ischr(c, exp      ) && possiblenum) 
      /* must be in a number field */
      + Bsign      *(ischr(c, sign     ) && (last_is & (Bexp | Bseparator | Beol)))
      /* must be after exponent or not in
      number because starting number */
      + Bcomment   * (ischr(c, comment  ) && strlen(comment))
      /* comment starts if we are waiting for it */
      + Bseparator * ischr(c, separator);
    
    
    if (is & Bseparator) c=' ';
    
    filestrcp[pos] = c;
    
    if (!(possiblecmt) && (is & Bcomment))  /* activate comment field */
    {
      possiblecmt = 1;
      startcmtpos = pos;
    }
    
    if (possiblecmt)
    {
      
      if (is & Beol)                        /* end comment on eol */
      {
        endcmtpos = pos;
        fieldend |= Bcomment;
        possiblecmt = 0;
        /*      if (endcharpos < startcharpos)
        endcharpos= startcmtpos - 1;
        if (startcharpos < endcharpos) fieldend |= Balpha;*/ 
      }
      else
      {
        is = last_is;
        filestrcp[pos] = ' ';
      }
    }
    
    if ( (!possiblecmt) && (!possiblenum) && (is & needforstartnum) && (last_is & (Bseparator | Beol)))    /* activate num search */
    {
      possiblenum = 1;
      startnumpos = pos;
      need = needforstartnum;
      inpoint = 0;
      inexp   = 0;
    }
    
    if (possiblenum && !(possiblecmt))     /* in num field */
    {
      found = is & need;
      
      /* last column : update when found and (EOL and not groupnum) or (EOL and groupnum and columns (not empty previous num line)) */
      /* OK : columns = 0 when EOL */
      /* OK : newline : when found after EOL (can be really before -> if columns == 1)*/
      /* end of num field : found && columns != last_columns */
      
      if ((last_is & (Bnumber | Bpoint)) && (is & (Beol | Bseparator)))
      {
        columns++; /* detects num end : one more column */
        if (found && (columns == 1)) rows++;  /* this is a new line starting */
      }
      
      if (is & Beol)
      {
        if (!groupnum)
          if ((columns != last_columns) && (startnumpos < last_eolpos))
          {                                    /* change in columns -> end of preceeding num field */
            endnumpos = last_eolpos - 1;
            pos       = last_eolpos;
            is        = Beol;
            need      = needforstartnum;
            fieldend |= Bnumber;
            endcharpos= startnumpos - 1;
            columns   = last_columns;
            if (startcharpos <= endcharpos) fieldend |= Balpha;
          }
          else
          {
            last_columns = columns;
            columns = 0;
          }
          else  /* group num */
          {
            if (!columns && last_columns)
            {
              filestrcp[pos] = ' ';  /* pass on : continue in same field */
            }
            else
            {
              if (columns && !last_columns)
              {
                last_columns = columns;  /* first line on matrix */
                columns = 0;
              }
              else
                if (last_columns && (columns != last_columns)  && (startnumpos < last_eolpos))
                {                                    /* change in columns -> end of preceeding num field */
                  endnumpos = last_eolpos - 1;
                  pos       = last_eolpos;
                  is        = Beol;
                  need      = needforstartnum;
                  fieldend |= Bnumber;
                  endcharpos= startnumpos - 1;
                  columns   = last_columns;
                  if (startcharpos <= endcharpos) fieldend |= Balpha;
                }
                else
                {
                  columns = 0;
                }
            }
          }
      }
      
      if (!found)
      {
        if (last_is & (Beol | Bseparator) && (need != Bnumber))  /* end of num field ok, except when the num started by a single point, not follwed by a 0-9 number */
        {
          endnumpos  = pos - 2;
          endcharpos = startnumpos - 1;
          if (last_columns == 0) last_columns = columns;
          if (startcharpos <= endcharpos) fieldend |= Balpha;
          fieldend  |= Bnumber;
        }
        else /* anomalous end of num */
        {
          if (startnumpos >= last_seppos) /* first possible number is not a number */
          {
            columns   = last_columns;
            possiblenum = 0; /* abort and pass */
            if (fieldend & Bnumber) fieldend -= Bnumber;
            possiblenum = 0;
          }
          else
            if ((columns > 0) && (startnumpos >= last_eolpos)) /* only a line */
            {
              endnumpos = last_seppos - 1;
              pos       = last_seppos;
              is        = Bseparator;
              need      = needforstartnum;
              fieldend |= Bnumber;
              endcharpos = startnumpos - 1;
              if (startcharpos <= endcharpos) fieldend |= Balpha;
            }
            else /* already passed more than one line */
            {
              endnumpos = last_eolpos - 1;
              pos       = last_eolpos;
              is        = Beol;
              need      = needforstartnum;
              fieldend |= Bnumber;
              endcharpos= startnumpos - 1;
              columns   = last_columns;
              if (startcharpos <= endcharpos) fieldend |= Balpha;
            }
        }
      }
      else                                /* still in num */
      {
        if (is & Bpoint)
        {
          if (inpoint || inexp) 
          {
            need = 0; 
          }
          else
          {
            if (last_is & (Bseparator | Beol)) need = Bnumber;
            else need = needafterpoint;
            inpoint = 1; 
            
          }
        }
        else
          if (is & Bsign)        need = needaftersign;
          else
            if (is & Bexp)       { need = needafterexp; inpoint = 0; inexp = 1; }
            else
              if (is & Bseparator) { need = needaftersep; inpoint = 0; inexp = 0; }
              else
                if (is & Bnumber)      need = needafternumber;
                else
                  if (is & Beol)       { need = needaftereol; inpoint = 0; inexp = 0; }
                  else
                    need = needafternumber;
      }
  }
  /* &&  (abs(endcharpos-startcharpos)+abs(endcharpos-startcharpos)>0) */
  if ((fieldend)   )
  {
    
    /* =========================== STRINGS/ASCII fields */
    if ((fieldend & Balpha) || (numheader && (fieldend & Bnumber)))
    {
      if (verbose >= 1)
      {
        printf("%5i %s %10i %10i %10i %10i : ",id,"chr",startcharpos,endcharpos,rows, columns);
        printfield(filestr, startcharpos, filesize, endcharpos - startcharpos);
        printf("\n");
      }
      if (tablextr)
      { 
        filetable[filetablecontent] = id;
        filetable[filetablecontent+1] = Balpha;
        filetable[filetablecontent+2] = startcharpos;
        filetable[filetablecontent+3] = endcharpos;
        filetable[filetablecontent+4] = rows;
        filetable[filetablecontent+5] = columns;
        filetablecontent += 6;
        if (filetablecontent > filetablesize)
        {
          filetablesize = filetablecontent+1024;
          filetable = (long *)realloc(filetable, 6*sizeof(long)*filetablesize);
          if (filetable == NULL)
          {
            printf("looktxt: Warn : unable to re-alloc %i bytes for filetable. No table.\n",6*sizeof(long)*filetablesize);
            fflush(stdout);
            tablextr = 0;
          }
        }
      }
      
      if ((!passascii || (numheader && (fieldend & Bnumber)))
        && (outtype == OUT_OCT) && (tablextr != 2) && (endcharpos - startcharpos + 1 > 0))
      {
        fprintf(fout,"# name: ");
        if (!noroot) fprintf(fout,"%s_",filename);
        fprintf(fout,"s%i\n",id);
        fprintf(fout,"# type: string\n");
        fprintf(fout,"# length: %i\n", endcharpos - startcharpos + 1);
        fwrite(filestr+startcharpos, abs(endcharpos - startcharpos + 1), 1, fout);
        fprintf(fout,"\n");
      }
      else
        if ((!passascii  || (numheader && (fieldend & Bnumber))) 
          && ((outtype == OUT_M))  && (tablextr != 2) && (endcharpos - startcharpos + 1 > 0))
        {
          if (outtype == OUT_M)
          {
            for (i=startcharpos; i <= endcharpos; i++)
            {
              if (filestr[i] == '\'')
                filestr[i] = '\"';
              else
              {
                if (!isprint(filestr[i]))
                {
                  if ((filestr[i] == '\n') || (filestr[i] == '\r'))
                    filestr[i] = '\t';
                  else
                    filestr[i] = ' '; 
                }
              }
            }  
            
            if (!noroot) fprintf(fout,"%s_",filename);
            if (usestruct) fprintf(fout,"s.");
            fprintf(fout,"s%i = '",id);
            fwrite(filestr+startcharpos, abs(endcharpos - startcharpos + 1), 1, fout);
            fprintf(fout,"';\n");
          }
        }
        
        startcharpos = pos;
        fieldend -= Balpha;
        if (!numheader) id++;
    }  
    /* && (abs(endnumpos-startnumpos)) */
    /* =========================== NUMERIC fields */
    if ((fieldend & Bnumber)  )
    {
      columns = last_columns;
      if (columns && (startnumpos <= endnumpos))
      {
        if (rows <= 0) rows = 1;
        if (verbose >= 1)
        {
          printf("%5i %s %10i %10i %10i %10i : ",id,"num",startnumpos,endnumpos,rows, columns);
          printfield(filestr, startnumpos, filesize, endnumpos - startnumpos);
          printf("\n");
        }
        if (tablextr)
        { 
          filetable[filetablecontent] = id;
          filetable[filetablecontent+1] = Bnumber;
          filetable[filetablecontent+2] = startnumpos;
          filetable[filetablecontent+3] = endnumpos;
          filetable[filetablecontent+4] = rows;
          filetable[filetablecontent+5] = columns;
          filetablecontent += 6;
          if (filetablecontent > filetablesize)
          {
            filetablesize = filetablecontent+1024;
            filetable = (long *)realloc(filetable, 6*sizeof(long)*filetablesize);
            if (filetable == NULL)
            {
              printf("looktxt: Warn : unable to re-alloc %i bytes for filetable. No table.\n",6*sizeof(long)*filetablesize);
              fflush(stdout);
              tablextr = 0;
            }
          }
        }
        
        prod = rows*columns;
        if (prod > maxprod) { maxprod = prod; maxfield = id; }
        if ((outtype == OUT_OCT) && (prod >= line)  && (tablextr != 2))
        {
          numextr++;
          fprintf(fout,"# name: ");
          if (!noroot) fprintf(fout,"%s_",filename);
          fprintf(fout,"n%i\n",id);
          if ((rows == 1) && (columns == 1))
            fprintf(fout,"# type: scalar\n");
          else
          {
            fprintf(fout,"# type: matrix\n");
            fprintf(fout,"# rows: %i\n", rows);
            fprintf(fout,"# columns: %i\n", columns);
          }
          for (i=startnumpos; i <= endnumpos; i++)
            if ( (i < startcmtpos) || (i > endcmtpos) )
              fprintf(fout,"%c", filestrcp[i]);
            fprintf(fout,"\n");
        }
        else
          if ((outtype == OUT_M) && (prod >= line)  && (tablextr != 2))
          {
            numextr++;
            if (!noroot) fprintf(fout,"%s_",filename);
            if (usestruct) fprintf(fout,"s.");
            fprintf(fout,"n%i = [ ",id);
            for (i=startnumpos; i <= endnumpos; i++)
              if ( (i < startcmtpos) || (i > endcmtpos) )
                fprintf(fout,"%c", filestrcp[i]);
              fprintf(fout," ];\n");
          } 
          else
            if ((outtype == OUT_TXT) && (prod >= line)  && (tablextr != 2))
            {
              numextr++;
              for (i=startnumpos; i <= endnumpos; i++)
                if ( (i < startcmtpos) || (i > endcmtpos) )
                  fprintf(fout,"%c", filestrcp[i]);
                fprintf(fout,"\n");
            }
      }
      possiblenum = 0;
      last_eolpos = endnumpos;
      startcharpos = endnumpos+1;
      fieldend -= Bnumber;
      columns = 0;
      rows = 0;
      last_columns = 0;
      if (i) pos = endnumpos+1;
      id++;
      inpoint = 0;
    }
    
    if (fieldend & Bcomment)
    {
      if (verbose >= 1)
      {
        printf("%5i %s %10i %10i %10i %10i : ",id,"cmt",startcmtpos,endcmtpos,rows, columns);
        printfield(filestr, startcmtpos, filesize, endcmtpos - startcmtpos);
        printf("\n");
      }
      if (tablextr)
      { 
        filetable[filetablecontent] = id;
        filetable[filetablecontent+1] = Bcomment;
        filetable[filetablecontent+2] = startcmtpos;
        filetable[filetablecontent+3] = endcmtpos;
        filetable[filetablecontent+4] = rows;
        filetable[filetablecontent+5] = columns;
        filetablecontent += 6;
        if (filetablecontent > filetablesize)
        {
          filetablesize = filetablecontent+1024;
          filetable = (long *)realloc(filetable, 6*sizeof(long)*filetablesize);
          if (filetable == NULL)
          {
            printf("looktxt: Warn : unable to re-alloc %i bytes for filetable. No table.\n",6*sizeof(double)*filetablesize);
            fflush(stdout);
            tablextr = 0;
          }
        }
      }
      
      possiblecmt = 0;
      fieldend -= Bcomment;
      startcharpos = pos+1;
      id++;
    }
  }
  
  if ((is & Beol))
  {
    last_eolpos = pos;
    last_seppos = pos;
  }
  
  if (is & Bseparator)
    last_seppos = pos;
  
  pos++;
}
while (pos <= filesize+1);

/* end scan */
if (verbose >= 0)
{
  printf("%s [%i fields, %i numerics] extracted into %s\n", filenamecp,id-1, numextr,fileout);
  if (numextr)
  {
    printf("Main numeric field is : ");
    if (!noroot && ((outtype == OUT_OCT) || (outtype == OUT_M))) printf("%s_",filename);
    if (usestruct && (outtype == OUT_M)) printf("s.");
    printf("n%i (%i elements)\n",maxfield, maxprod);
  }
}
fflush(NULL);
free(filestr); 
free(fileout);
free(filestrcp);
free(filenamecp);

if (tablextr)
{
  if ((outtype == OUT_M) || (outtype == OUT_OCT))
  {
    fprintf(fout,"%% Here follows the table matrix\n");
    fprintf(fout,"%% It describes how was analyzed the file %s\n", filename);
    fprintf(fout,"%% and shows all blocks found in it.\n");
    fprintf(fout,"%% [ id type start_pos end_pos rows columns ]\n");
  }
  if (outtype == OUT_M)
  {
    if (!noroot) fprintf(fout,"%s_",filename);
    if (usestruct) fprintf(fout,"s.");
    fprintf(fout,"table = [ ");
  }
  if (outtype == OUT_OCT)
  {  
    fprintf(fout,"# name: ");
    if (!noroot) fprintf(fout,"%s_",filename);
    fprintf(fout,"table\n");
    fprintf(fout,"# type: matrix\n");
    fprintf(fout,"# rows: %i\n", id-2);
    fprintf(fout,"# columns: %i\n", 6);
  }
  for (i = 0; i < filetablecontent; i += 6)
    fprintf(fout,"%5i %2i %10i %10i %10i %10i\n",filetable[i],filetable[i+1],filetable[i+2],filetable[i+3], filetable[i+4], filetable[i+5]);
  
  if (outtype == OUT_M) fprintf(fout," ];\n");
  
  free(filetable);
  if (verbose >= 1)
  {
    printf("Table of fields is : ");
    if (!noroot) printf("%s_",filename);
    if (usestruct && (outtype == OUT_M)) printf("s.");
    printf("table\n");
  }
}
if ((outtype == OUT_M) || (outtype == OUT_OCT)) { fprintf(fout,"%% End of file %s\n",filename); fclose(fout); }
return (maxfield);
}

/* MAIN ====================================================== */

long main(int argc, char *argv[])
{
  long  i;
  char *arg;
  char *tmp;
  
  verbose  = 0;
  passascii= 1;
  groupnum = 0;
  noroot   = 0;
  outtype  = OUT_M;
  forcefile= 0;
  line     = 0;
  tablextr = 0;
  username = 0;
  
  
  if (argc == 1)
  {
    printf("Usage : looktxt [options] filename\n");
    printf("Action: Search and export numerics in a text/ascii file.\n");
    printf("Type looktxt -h for extended help, looktxt -v to get default parameters.\n");
    printf("Example : looktxt foo.txt \n");
    printf(author " version " version " (" date ")\n");
    return(-1);
  }
  else
  {     
    arg = (char *)malloc(1024);
    point = (char *)malloc(256);
    separator = (char *)malloc(256);
    eol = (char *)malloc(256);
    comment = (char *)malloc(256);
    
    strcpy(point,Cpoint);
    strcpy(separator,Cseparator);
    strcpy(eol,Ceol);
    strcpy(comment,Ccomment);
    
    userfile = (char *)malloc(1024);
    filename = (char *)malloc(1024);
    filename[0] = '\0';
    userfile[0] = '\0';
    
    for (i=1; i < argc; i++)
    {
      strcpy(arg,argv[i]);
      
      if (strlen(arg) > 255)
      {
        printf("looktxt: arg %i too long\n",i-1);
        return(-1);
      }
      if (arg[0] == '-')  /* option ? */
      {
        if (strncmp(arg,"-v",2) == 0)
          verbose = 1;
        else
          if (strncmp(arg,"-D",2) == 0)
            verbose = 2;
          else
            if (strncmp(arg,"-H",2) == 0)
              numheader = 1;
            else
              if (strncmp(arg,"-S",2) == 0)
                usestruct = 1;
              else
                if (strncmp(arg,"-p",2) == 0)
                {
                  if (strlen(arg)>=4)
                    strcpy(point,arg+3);
                  else printf("looktxt: warning: no point ! using default.\n"); 
                }
                else
                  if (strncmp(arg,"--silent",8) == 0)
                  {
                    verbose = -1;
                  }
                  else
                    if (strncmp(arg,"-c",2) == 0)
                    {
                      if (strlen(arg)>=4)
                        strcpy(comment,arg+3);
                      else strcpy(comment,"");
                    }
                    else
                      if (strncmp(arg,"-s",2) == 0)
                      {    
                        if (strlen(arg) >= 4)
                          tmp = (char *)make_esc_char((char *)(arg + 3));  
                        else tmp = NULL;
                        if (tmp != NULL)
                        { free(separator); separator = tmp; }
                        else printf("looktxt: warning: no separator ! using default.\n");
                      }
                      else
                        if (strncmp(arg,"-n",2) == 0)
                        {
                          if (strlen(arg) >= 4)
                            tmp = (char *)make_esc_char((char *)(arg+3));
                          else tmp = NULL;
                          if (tmp != NULL)
                          { free(eol); eol = tmp; }
                          else printf("looktxt: warning: no end of line ! using default.\n");
                        }
                        else
                          if (strncmp(arg,"-F",2) == 0)
                          {
                            if (strlen(arg) >= 4)
                            { username = 1;
                            strcpy(userfile,arg+3); }
                            else printf("looktxt: warning: no outputfile ! using default.\n");
                          }
                          else
                            if ((strncmp(arg,"-o",2) == 0) & (strlen(arg) > 3))
                            {
                              if (arg[3] == 'o')
                                outtype = OUT_OCT;
                              else
                                if (arg[3] == 'm')
                                  outtype = OUT_M;
                                else
                                  if (arg[3] == 't')
                                    outtype = OUT_TXT;
                                  else
                                    outtype = OUT_M;
                            }
                            else
                              if (strncmp(arg,"-a",2) == 0)
                                passascii = 0;
                              else
                                if (strncmp(arg,"-g",2) == 0)
                                  groupnum = 1;
                                else
                                  if (strncmp(arg,"-f",2) == 0)
                                    forcefile = 1;
                                  else
                                    if (strncmp(arg,"-r",2) == 0)
                                      noroot = 1;
                                    else
                                      if ((strncmp(arg,"-l",2) == 0) & (strlen(arg) > 3))  /* -l=X */
                                      {
                                        line = atoi(arg+3);
                                        if (line < 0) line = 0;
                                        if (line > 60000) line = 60000; 
                                      }
                                      else
                                        if (strncmp(arg,"-t",2) == 0)
                                          tablextr = 1;
                                        else
                                          if (strncmp(arg,"-T",2) == 0)
                                            tablextr = 2;
                                          else
                                            if (strncmp(arg,"-h",2) == 0)
                                            {
                                              printf("Usage : looktxt [options] filename\n\n");
                                              printf("Action: Search and export numerics in a text/ascii file.\n");
                                              printf("   The programs looks into your file some numeric fields (and optionally\n");
                                              printf("   characters). Each identified numeric field is then named and exported\n");
                                              printf("   into an output filename. The output file can be a .m (Matlab), .oct\n");
                                              printf("   (Octave) or .txt file, with automatic or specified user name. For a\n");
                                              printf("   special type of data, you can specify point, separator, end of line\n");
                                              printf("   and comment conventions. A 'table' field can be exported, containing\n");
                                              printf("   informations about all other found fields (see end of this help).\n");
                                              printf("Examples : \n");
                                              printf("   'looktxt foo.txt'\n");
                                              printf("      will create 'foo_txt.m' file with all numerics found in 'foo.txt'.\n");
                                              printf("      Fields have names such as 'foo_txt_n12'\n");
                                              printf("   'looktxt -f -F=\"mynums\" -l=100 -v foo.txt' \n");
                                              printf("      will create/overwrite the file 'mynums.m' with only\n");
                                              printf("      vectors/matrices containing at least 100 elements, and showing all\n");
                                              printf("      the process and parameters.\n");
                                              printf("   'looktxt -a -g -r -o=\"oct\" foo.txt'\n");
                                              printf("      will search character and numeric fields, grouping them if\n");
                                              printf("      possible, naming them with short names (such as 'n12' or 's10'),\n");
                                              printf("      and creating 'foo_txt.oct' file\n");
                                              printf("Options :\n");
                                              printf("  [ -a ]               to look also for char fields (\\n and \\r -> \\v).\n");
                                              printf("  [ -c=\"{#|%%|c}\" ]     to set comment start character (till EOF)\n");
                                              printf("  [ -F=\"filename\" ]    to force output filename to 'filename.extension'\n");
                                              printf("  [ -f ]               to force creation/overwriting of output file\n");
                                              printf("  [ -g ]               to group numeric fields when possible\n");
                                              printf("  [ -H ]               extract ascii header preceeding numerics\n");
                                              printf("  [ -h ]               this help\n");
                                              printf("  [ -l=X ]             to extract only num fields with at least X elements\n");
                                              printf("  [ -n=\"{\\n|\\f}\" ]     to set end of line equivalent character\n");
                                              printf("  [ -o=\"{oct|m|txt}\" ] to set output type (ext) .oct, .m (default) or .txt\n");
                                              printf("  [ -p=\"{.|,}\" ]       to set point equivalents\n");
                                              printf("  [ -r ]               to use short field names \n");
                                              printf("                       (such as 'nXX' instead of 'filename_nXX')\n");
                                              printf("  [ -S ]               extract as a Matlab structure (with -o=\"m\")\n");
                                              printf("  [ -s=\"{\\t|\\v|\\r|,| |;}\" ] to set numeric separators\n");
                                              printf("  [ -t ]               to get field info 'table' matrix\n");
                                              printf("  [ -T ]               to extract only 'table' matrix\n");
                                              printf("  [ -v ]               to set verbose mode (show process and parameters)\n");
                                              printf("  [ --silent ]         not to print any message except errors\n");
                                              printf("The 'table' field (with -t option) is a matrix of rows :\n");
                                              printf("      [ id type start_pos end_pos rows columns ]\n");
                                              printf("describing all exported fields (identificator, type, start, end,\n");
                                              printf("dimensions) with type = 1 (num) | 2 (text) | 64 (comment)\n\n");
                                              printf("Type : looktxt -v to get default parameters.\n");
                                              printf(author " version " version " (" date ")\n");
                                              if (verbose == 2) { printf("main:free in help point separator eol comment filename userfile arg\n"); fflush(stdout); }
                                              free(point);
                                              free(separator);
                                              free(eol);
                                              free(comment);
                                              
#ifndef MATLAB_MEX_FILE
                                              free(filename);
                                              free(userfile); 
#endif
                                              
                                              free(arg);
                                              return(-1);
                                            }
                                            else
                                            {
                                              printf("looktxt: unknown option %s\n",arg);
                                            }
      } /* if option */
      else
      { 
        strcpy(filename,arg);
      }
    }
    if (verbose >= 1)
    {
      printf("* looktxt " version " (" date ")");
      if (verbose == 1)
        printf(": verbose mode\n");
      else
        printf(": DEBUG mode\n");
      printf("filename = \"%s\"\n",filename);
      printf("-p point = \"%s\"\n",point);
      printf("-s separator = \"%s\"\n",undo_esc_char(separator));
      printf("-n eol = \"%s\"\n",undo_esc_char(eol));
      printf("-c comment = \"%s\"\n",comment);
      printf("-l Extract num fields with more than %i elements\n",line);
      if (passascii) printf("-a Pass ascii mode\n");
      if (groupnum) printf("-g Group numerics mode\n");
      if (noroot) printf("-r No root for field names mode\n");
      if (outtype == OUT_OCT) printf("-o Output is .oct file\n");
      if (outtype == OUT_M) printf("-o Output is .m file\n");
      if (outtype == OUT_TXT) printf("-o Output is .txt file\n");
      if (forcefile) printf("-f Force creation of output file\n");
      if (tablextr == 1) printf("-t Extract field info in table 'matrix'\n");
      if (usestruct) printf("-S Use Matlab structure\n");
      if (numheader) printf("-H Extract ascii header before each num\n");
      if (username) printf("-F User output root filename \"%s\"\n",userfile);
    }
    if (strlen(filename) > 0)
    {
      i=scanfile(filename);
    }
    else
    { printf("looktxt: no file to process (looktxt -h for help)\n"); i = -2; }
    
    free(point);
    free(separator);
    free(eol);
    free(comment);
#ifndef MATLAB_MEX_FILE
    free(filename);
    free(userfile); 
#endif
    free(arg);
  } /* for arg */
  return(i);
}
