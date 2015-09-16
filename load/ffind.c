/********************************************************************* 
 ffind.c 	Usage fpos=ffind(filename,string)

 This mex file is intended to help matlab deal with large data files.
 It returns a pointer on the 'string' token in file 'filename'.
 An optional third parameter enables scanning of all 'string' positions.
 The matlab function fseek can be used.
 
 M. Zinkin 9.5.94
**********************************************************************/

#include <stdlib.h>
#include <stdio.h> 
#include <string.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include "mex.h"
                                                                  
#define MAXLINE   1024   
#define READFILE  prhs[0]
#define STRING    prhs[1] 
#define OUT       plhs[0]

/* version 5 -> 4 modifs :  mxArray -> Matrix, mxREAL -> 0 */

void   mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
                    
void   mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{                                          
  FILE *ReadPtr;
  char ReadFile[MAXLINE+1], String[MAXLINE+1];        
  double fpos;

  double *PosPtr;
  long PosNb, PosBuf;
  struct stat stfile;
  long filesize;
  char *filestr, *Line;
  char extractall;
  char *foundok;

  PosNb  = 0;
  PosBuf = 1024;

  if (nrhs!=2 && nrhs!=3) 
  {
    mexPrintf("Syntax:  <filepos>=ffind(<filename>,<search string>,{extract all flag}\n)");
    mexPrintf("         The last input parameter is optional\n");
    mexErrMsgTxt("Example: pos=ffind('myfile.txt','search','extract all');");
  }

  if (nlhs > 1)
    mexErrMsgTxt("ffind: Too many output arguments.");
    
  /* Check data type of input argument */
  if (!(mxIsChar(prhs[0]))) 
    mexErrMsgTxt("ffind: Input 1 must be of type char (string).");

  if (!(mxIsChar(prhs[1]))) 
    mexErrMsgTxt("ffind: Input 2 must be of type char (string).");

  if (nrhs == 3) 
    extractall = 1;
  else
    extractall = 0;

	mxGetString(READFILE,ReadFile,MAXLINE);
	mxGetString(STRING,String,MAXLINE);
  if (strlen(String) == 0)
    mexErrMsgTxt("ffind: Input 2 must be a non empty string.");


	
	PosPtr = (double*)mxMalloc(PosBuf*sizeof(double));

  if(stat(ReadFile,&stfile) != 0)
    mexErrMsgTxt("ffind: Fatal: Data file not found (does not exist).");
  filesize = stfile.st_size;
  filestr = (char *)mxMalloc(filesize+1);
  if (filestr == NULL)
    mexErrMsgTxt("ffind: Fatal : Unable to malloc for filestr (no memory left)");

/* Open file to read */
  if ((ReadPtr = fopen(ReadFile,"r"))==NULL) 
    mexErrMsgTxt("ffind: Fatal: Data file can not be opened (permissions ?).");              

/* read all */
  if (!fread(filestr, 1, filesize, ReadPtr))
    mexErrMsgTxt("ffind: Fatal : Unable to read file (permissions ?)");
  fclose(ReadPtr);

  filestr[filesize] = '\0';
  Line = filestr;
  fpos = -1;
  foundok = NULL;

	/* Read through file until find Start */
	do 
  {
    foundok = strstr(Line,String); /* !NULL when found */
    if (foundok != NULL)
    {
      Line = foundok+strlen(String);
      fpos = (double)(foundok - filestr);
      PosPtr[PosNb] = (double)fpos;
      PosNb++;
      if (PosNb >= PosBuf)
      {
        PosBuf += 1024;
        PosPtr = mxRealloc(PosPtr, PosBuf*sizeof(double));
      }
    } 
  }
	while ((foundok != NULL) && (extractall)); 
/* while ((fgets(Line,MAXLINE,ReadPtr)!=NULL) && (strstr(Line,String)==NULL)); */
	   
	if (fpos == -1) PosPtr[PosNb] = (double)fpos;
  mxFree(filestr);

	/* Return position in file */
	OUT=mxCreateDoubleMatrix(1,PosNb,0);
  PosPtr = mxRealloc(PosPtr, PosNb*sizeof(double));
	mxSetPr(OUT, PosPtr);
}
