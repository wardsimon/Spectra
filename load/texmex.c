/*-----------------------------program texmex.c-----------------------------*/

/*
% [mainfield, rootname, output, filestr] = looktxt('filename [options]') Import text data
% Converts any text file into a Matlab text file
% Then, type 'filename' to import data.
% 'rootname' is the root name of fields extracted to output
% 'output'   is the output file name (extension .m is to be added).
% 'filestr'  is the output filename contents for direct evaluation
%  
% Use "looktxt" or "looktxt -h" for more help about options. 
% Example :
%          looktxt foo.txt -p="." -s="\t\v,;" -c="#"
% E.Farhi 01/97 and K. Dabertran 04/99. v0.31 (05/04/01)
*/

/* 
% In the mexFunction , outputs arguments are:
% * the first variable corresponds to the main numeric field number.
% * the 2nd represents the stem of the variable name.
% * the 3rd represents the stem of the file name.
% * the 4th is the string of the output file for direct evaluation

% Author:  EF <farhi@ill.fr>
% Description: Import any text file data. Need 'looktxt.c' source code.

 * compile with : mex -O -output looktxt texmex.c
 *           or : mex -v -argcheck -output looktxt texmex.c 
 * content: C language, MEX library functions
 * tab = 2 chars
 */

#include <mex.h>	/* include MEX library for Matlab */

#define printf mexPrintf	/* Addapt looktxt.c code to Mex syntax */
#define malloc mxMalloc
#define realloc mxRealloc
#define version "0.9.4 (MeX 0.31)"
/* #define free mxFree  */
#define free NoOp 

#define main looktxt		/* change stand-alone source to a Matlab usable library */
#define argc carg
#define argv varg
#define TEXMEX

#include "looktxt.c" 		/* makes all the job */

int NoOp(char *pointer)
{
  return 0;
}

void mexFunction(int nlhs, mxArray *plhs[], 
                        int nrhs, const mxArray *prhs[])
{
  char *InputTokens, *InputString;
  char *varg[128];
  char EndFlag = 0;
  int  carg, i, mrows, ncols;
  int  status;
  int  buflen;
  double *pMainField;
  long MainField;
  char lexeme[1024];
  char *StartLexeme, *EndLexeme;
  char *EndString;
  
  struct stat stfile;
  FILE *fout;
  char *filestr;
  long filesize;
  
  passascii= 1;
  groupnum = 0;
  noroot   = 0;
  outtype  = OUT_M;
  forcefile= 0;
  line     = 0;
  tablextr = 0;
  username = 0;
  numheader= 0;
  usestruct= 0;
  verbose  = 0;

  carg    = 1;
  buflen = 1024;
  varg[0] = (char*)mxMalloc(150);
  strcpy(varg[0],"looktxt");

  /* check in/out parameters */

  if (nlhs > 4)
	  mexErrMsgTxt("looktxt : Too many output arguments (4 max).");

  /* allocate memory */
  
  for (i = 0; i < nrhs; i++)
  {
    if (mxIsChar(prhs[i]) != 1)
    {
      mexPrintf("looktxt/mex : argument %i\n", i);
      mexErrMsgTxt("looktxt/mex : Input should be strings");
    }
      
    buflen      = (mxGetM(prhs[i])*mxGetN(prhs[i]))+1;
    InputString = (char*)mxMalloc(buflen+64);
    if (InputString == NULL)
    {
      mexPrintf("looktxt/mex : argument %i. Size %i\n", i, buflen);
      mexErrMsgTxt("looktxt/mex : can not allocate memory for input string\n");
    }
    status      = mxGetString(prhs[i], InputString, buflen);
    if (status != 0)
    {
      mexPrintf("looktxt/mex : argument %i. Status %i\n", i, status);
      mexErrMsgTxt("looktxt/mex : can not get input parameter\n");
    }

    /* cut input string into separated arguments for main(argc,argv) syntax */ 

    EndFlag = 0;
    StartLexeme = InputString;
    EndString = InputString+strlen(InputString);
    
    while ((EndFlag == 0) && (carg < 128))
    {
      if (*StartLexeme == ' ')
        while (*StartLexeme == ' ' && StartLexeme < EndString)
        {       /* look for first non ' ' : StartLexeme */
          StartLexeme++;  /* pass all spaces */
        }
      /* now StartLexeme points on a non space or is at end */
        
      if (*StartLexeme == '\0' || StartLexeme >= EndString) EndFlag = 1;
      else
      {
        /* look for position of first next ' ' : EndLexeme */
        EndLexeme = strchr(StartLexeme+1, ' ');

        if (EndLexeme == NULL)
          EndLexeme = EndString;
          
        if (EndLexeme - StartLexeme > 0)
        {
          strncpy(lexeme, StartLexeme, EndLexeme - StartLexeme+1);
          lexeme[EndLexeme - StartLexeme] = '\0';
        
          /* lexeme = (char *)strtok(InputTokens, " ,");
          InputTokens = NULL; */
          /* replacement to strtok that has internal non protected malloc ! */
          StartLexeme = EndLexeme+1;
        }
      }
      
      if (strlen(lexeme) != 0 && lexeme != NULL && EndFlag == 0)
      {
        varg[carg] = (char*)mxMalloc(strlen(lexeme)+64);
        strcpy(varg[carg], lexeme);
        carg++;
      }
      else 
        EndFlag = 1;
    }
    
  } /* end for nrhs */
  
/*  for (i=0; i < carg; i++)
    mexPrintf("looktxt/mex: arg %i = <%s>\n", i, varg[i]); */

  /* link integer to matlab double */

  MainField = looktxt(carg,varg);  /* call main routine */
  

    /* set output parameters */
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    pMainField = mxGetPr(plhs[0]);
    *pMainField = (double)MainField;
    
  if (MainField != -1)
  {

    if (strlen(filename))
      plhs[1] = mxCreateString(filename);
    if (strlen(userfile))
      plhs[2] = mxCreateString(userfile);

    /* now opens output 'userfile' file if required, to export it into matlab string array */
    if ((nlhs > 3) && (strlen(userfile)))
    {
      EndFlag = 0;
      userfile = (char *)strcat(userfile,".m");
      fout = fopen(userfile, "r");
      if (fout == NULL)
      {
        mexPrintf("looktxt/mex: Warning : unable to open file %s\n", userfile); 
        EndFlag = 1; 
      }
      if (EndFlag == 0)
      {
        stat(userfile, &stfile);
        filesize = stfile.st_size;
        filestr  = (char*)mxMalloc(filesize+64);
        if (filestr == NULL)
        {
          mexPrintf("looktxt/mex: Warning : unable to malloc %i bytes for filestr (%s)\n", filesize, userfile);
          EndFlag = 2;
        }
      }
      if (EndFlag == 0)
      {
        if (!fread(filestr, 1, filesize, fout))
        {
          mexPrintf("looktxt/mex: Warning : unable to read file %s\n", userfile);
          EndFlag = 3;
        }
      }
      if (EndFlag != 1)
        fclose(fout);
      if ((EndFlag != 2) && (EndFlag != 1))
      {
        filestr[filesize] = '\n';
        filestr[filesize+1] = '\0';
        plhs[3] = mxCreateString(filestr);
        free(filestr);
      }
      else
        plhs[3] = mxCreateString("");
    }
  }
  

  for (i=0; i < carg; free(varg[i++]));
}
