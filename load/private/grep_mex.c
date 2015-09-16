/*
 * grep_mex.c       - Regular expression matching
 * 
 * C-mexfunction for regular expression matching.
 * [i_start, i_end] = grep_mex(s, pattern, options)
 * s: string to search in
 * pattern: regular expression pattern
 * options: optional parameter that may contain
 *  'i' for 'ignore case
 *  'o' to only get the first match
 * 
 * Ref: man regex
 *      http://www.gnu.org/manual
 * 
 * See also: Grep
 *
 * Part of: utilities load/private
 * Author:  E. Farhi <farhi@ill.fr>. Feb 6th, 2003.
 *
 * Used in: Grep
 * Based on 'brede/bredx_re_match.c' from Lyngby and Brede Matlab softwares
 * at <http://hendrix.imm.dtu.dk/software/lyngby>
 */

#include <sys/types.h>
#include <regex.h>
#include <mex.h> 

/* Input arguments */
#define S        prhs[0]
#define P        prhs[1]
#define Options  prhs[2]

/* Output arguments */
#define M        plhs[0]
#define N        plhs[1]


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  
  double  *matchs;                          /* The Result */
  unsigned int    dS, dP, dOptions;


  /* regex associated variables */
#define ERRBUF 128
#define PREG 1024
  regex_t preg[PREG];      /* Pointer to pattern buffer */
  char *regex;             /* Pointer to the null-terminated string */
  char *string;            /* String to be searched */
  char *options;           /* Options for search (i,o) */
  int cflags;              /* Flags used to determine the type of compilation */
  int errcode;             /* Return code from regcomp */
  int retcode;             /* Return code from regexec */
  char errbuf[ERRBUF];    

  size_t nmatch;
  regmatch_t pmatch;
  int eflags;
  int pindex=0;
  int offset=0;

  int buflenP, buflenS, buflenO;
  int start[PREG];
  int end[PREG];

  /* Check for proper number of arguments */ 
  if (nrhs < 2) mexErrMsgTxt("grep_mex: Too few input arguments (min 2)."); 
  if (nrhs > 3) mexErrMsgTxt("grep_mex: Too many input arguments (max 3)."); 
  if (nlhs > 2) mexErrMsgTxt("grep_mex: Too many output arguments (max 2)."); 

  /* Check the dimensions of Input par. */
  dS       = mxGetNumberOfDimensions(S);
  dP       = mxGetNumberOfDimensions(P);
  if (dS > 2) mexErrMsgTxt("grep_mex: String must be a single char line.");
  if (dP > 2) mexErrMsgTxt("grep_mex: RegExp must be a single char line.");
  if (!mxIsChar(S)) mexErrMsgTxt("grep_mex: S should be a char array."); 
  if (!mxIsChar(P)) mexErrMsgTxt("grep_mex: P should be a char array."); 

  cflags = REG_EXTENDED | REG_NEWLINE;
  if (nrhs >= 3) {
    if (!mxIsChar(Options)) mexErrMsgTxt("grep_mex: Options should be a char array."); 
    dOptions = mxGetNumberOfDimensions(Options);
    if (dOptions > 2) mexErrMsgTxt("grep_mex: Options must be a single char line.");
    buflenO = (mxGetM(Options) * mxGetN(Options) * sizeof(mxChar)) + 1;
    options = malloc(buflenO);
    if (!options) mexErrMsgTxt("grep_mex: can not allocate Options."); 
    mxGetString(Options, options, buflenO);
    if (strchr(options, 'i')) cflags |= REG_ICASE;
    if (strchr(options, 'o')) pindex = 1;
  }

  /*  Copy a string mxArray's data into a C-style string */
  buflenS = (mxGetM(S) * mxGetN(S) * sizeof(mxChar)) + 1;
  buflenP = (mxGetM(P) * mxGetN(P) * sizeof(mxChar)) + 1;
  string = malloc(buflenS);
  regex = malloc(buflenP);
  mxGetString(S, string, buflenS);
  mxGetString(P, regex, buflenP);
  
  /* compile the pattern */
  if (errcode = regcomp(preg, regex, cflags)) {
    regerror(errcode, preg, errbuf, ERRBUF);
    mexErrMsgTxt(errbuf);     
  }
  
  /* main loop for searching pattern */
  offset=0;
  nmatch=0;
  do {
    /* Do the pattern matching */
    retcode = regexec(preg, string+offset, 1,  &pmatch, eflags);
    if (retcode != REG_NOMATCH)
    { /* store start and end indexes */
      if (pmatch.rm_so >= 0 && pmatch.rm_eo >= pmatch.rm_so) {
        start[nmatch] = offset+pmatch.rm_so+1;
        end[nmatch]   = offset+pmatch.rm_eo;
        nmatch++;
        if (pmatch.rm_eo) offset+=pmatch.rm_eo; else offset++;
      } else retcode = REG_NOMATCH;
    }
  } while (!retcode && offset < strlen(string) && nmatch < PREG && !pindex);

  /* Free the memory allocated */
  regfree(preg);
  free(regex);
  free(string);

  /* regexec returns zero for a successful match or REG_NOMATCH for failure. */
  if (retcode && retcode != REG_NOMATCH) {
    mexErrMsgTxt("Internal error: Could not understand return argument from regexec");
  }

  /* Create a matrix for the return argument */  
  M  = mxCreateDoubleMatrix(1, nmatch, mxREAL);
  N  = mxCreateDoubleMatrix(1, nmatch, mxREAL);
  
  /* Fill up return argument */
  for (pindex=0; pindex < nmatch; pindex++)
  {
    mxGetPr(M)[pindex] = (double) (start[pindex]);
    mxGetPr(N)[pindex] = (double) (end[pindex]);
  }

  return;
}

