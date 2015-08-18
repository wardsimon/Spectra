#include <math.h>
#include <stdio.h>
#include "mex.h"

/* Input Arguments */

#define	NMC_IN	prhs[0]


/* Output Arguments */

#define	S_OUT	plhs[0]

#define pi 3.14159265  

float gasdev(long *idum, double r[], int nmc) 
{
	float ran1(long *idum);
	static int iset=0,i;
	static float gset;
	float fac,rsq,v1,v2;

    for (i=0; i<nmc; i++)
    {
        if  (iset == 0) {
            do {
                v1=2.0*ran1(idum)-1.0;
                v2=2.0*ran1(idum)-1.0;
                rsq=v1*v1+v2*v2;
            } while (rsq >= 1.0 || rsq == 0.0);
            fac=sqrt(-2.0*log(rsq)/rsq);
            gset=v1*fac;
            iset=1;
            r[i] = v2*fac;
        } else {
            iset=0;
            r[i]=gset;
        }
    }
    return;
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum) 
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


#ifdef __STDC__
void mexFunction(
	int		nlhs,
	mxArray	*plhs[],
	int		nrhs,
	const mxArray	*prhs[]
	)
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[];
const mxArray *prhs[];
#endif  
{
	double	*r;
	double	*nmc;
    long idum=10634;
	if (nrhs != 1) {
		mexErrMsgTxt("ran1 requires one input arguments.");
	}

	/* Assign pointers to the various parameters */
	nmc = mxGetPr(NMC_IN);

    /* Create a matrix for the return argument */

	S_OUT = mxCreateDoubleMatrix(1, (int)nmc[0], mxREAL);
    r = mxGetPr(S_OUT);

	/* Do the actual computations in a subroutine */

	gasdev(&idum,r,(int)nmc[0]);
	return;
}



                                          

