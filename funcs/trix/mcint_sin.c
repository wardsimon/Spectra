  /*

  TRIXMEX.C	

	The calling syntax is:

  	[s]=mcint(nmc,hkle,d,fwhm,b)

  	D.A. Tennant 29th February 1996

	compilation: cc mcint.c /g_float/include=matlab:[extern.include]
        on av6       cmexg mcint

	On PC platforms, must be compiled with /Alfw and /Gs option to generate
	proper code.

*/

/* to be done: pass variables in main over as mex and check defn's
 of hkle &c. */

#include <math.h>
#include <stdio.h>
#include "mex.h"

/* Input Arguments */

#define	NMC_IN	prhs[0]
#define	HKLE_IN	prhs[1]
#define D_IN	prhs[2]
#define FWHM_IN	prhs[3]
#define B_IN	prhs[4]


/* Output Arguments */

#define	S_OUT	plhs[0]

#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))

/* Define Constants */

#define pi 3.14159265     
#define bf 11.609        /* Bose factor conversion kT->meV */
#define fwsig 1 /* = 2.35482005 Conversion from FWHM to sigma =2*sqrt(2*ln(2)) */

/* Define the Monte Carlo routine. */
static
/* #ifdef __STDC__ */
void mcint(
	double s[],
	double nmc[],
	double hkle[],
	double d[],
	double fwhm[],
	double b[]
	)
/* #else
mcint(s,nmc,hkle,d,fwhm,b)
double s[];
double nmc[];
double hkle[],d[],fwhm[],double b[];
#endif */
{
	double cs(double hkle[], double d[]); /* called func. */
	float gasdev(long *idum);                             /* called func. */
	long idum=10634;         /* seed long integer for random no. gen. */
	float ans;
	int i,j,index;
	
        
/* fwhm of principle axes and transform matrix */

	double xx[4],q[4],b44[4][4],ans1,r1,r2,r3;
	double yinc=0,y=0,sumcs2=0;
	index=0;

	for (i=0;i<4; ++i){
		for (j=0;j<4; ++j){
			b44[i][j]=b[index++];
		}
 	}                                   
	for (i=0; i<nmc[0]; ++i) {  
        	for (j=0;j<4; ++j){
			r1=gasdev(&idum);
			xx[j]=r1*fwhm[j]/fwsig;
	        } 
		q[0]=hkle[0];
		q[1]=hkle[1];
		q[2]=hkle[2];
		q[3]=hkle[3];
		for (j=0;j<4; ++j){
			q[0] += b44[0][j]*xx[j]; 
			q[1] += b44[1][j]*xx[j]; 
			q[2] += b44[2][j]*xx[j]; 
			q[3] += b44[3][j]*xx[j]; 
	        } 

		yinc=cs(q,d);
		y += yinc;
		sumcs2 += yinc*yinc; 
/*        	printf("dispersion is %6.3f \n",ans1);       */
        }
/* Check the uncertainty in the integral and do more Monte Carlo points	
   if necessary. */
        s[0]=y/nmc[0];
	return;  
}




/* Now define the cross-section and dispersion */

/* to be made into a separate linkable routine */
/* parameters for the cross section calc.      */
/* disp = dispersion parameters                */
/* d[0]  = width	            */
/* d[1]  = temperature	        */
/* d[2]  = BZ centre for H      */
/* d[3]  = BZ centre for K      */
/* d[4]  = BZ centre for L      */
/* d[5]  = Periodicity along H  */
/* d[6]  = Periodicity along K  */
/* d[7]  = Periodicity along L  */
/* d[8]  = BZ boundary energy H */
/* d[9]  = BZ boundary energy K */
/* d[10] = BZ boundary energy L */
/* d[11] = Gap energy           */

double cs(
	double  hkle[],
	double	d[]
	)
{
	double	wqs,qx,qy,qz,Ah,Ak,Al;
	double  nw,en,s;
	/* dispersion calculation first! */
	/* position in reduced q */
	qx=hkle[0]-d[2]; 
	qy=hkle[1]-d[3];
	qz=hkle[2]-d[4];

	/* Amplitudes along q */
	Ah=d[8]*d[8]-d[11]*d[11];
	Ak=d[9]*d[9]-d[11]*d[11];
	Al=d[10]*d[10]-d[11]*d[11];
	
/*  Account for the gap   */
	en = hkle[3]+d[11];
/*	Now the dispersion, just a simple sin function*/
	wqs = Ah*sin(d[5]*pi*qx)*sin(d[5]*pi*qx);  
	wqs = wqs + Ak*sin(d[6]*pi*qy)*sin(d[6]*pi*qy);  
	wqs = wqs + Al*sin(d[7]*pi*qz)*sin(d[7]*pi*qz);
	wqs = wqs + d[11]*d[11];  
	
/*	mexPrintf("%f %f %f %f  \n",qx,qy,qz,wqs);	*/	

        /* now cross-section part! - Damped Harmonic Oscillator. */

	nw=(4.0/pi)*d[0]*hkle[3]/(1-exp(-hkle[3]*bf/d[1])); /* temperature factor. */
	s=nw/((hkle[3]*hkle[3]-wqs)*(hkle[3]*hkle[3]-wqs)+4*hkle[3]*hkle[3]*d[0]*d[0]); 
	return s;
}   
/* end of cross section routine */
/* ____________________________ */


/* Numerical recipies random number generators */
/* #include <math.h> */

float gasdev(long *idum) 
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software #17Z9'2150. */

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
/* (C) Copr. 1986-92 Numerical Recipes Software #17Z9'2150. */

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
	double	*s;
	double	*nmc;
	double	*hkle,*d,*fwhm,*b;
	unsigned int	m,n;

	/* Check for proper number of arguments */

	if (nrhs != 5) {
		mexErrMsgTxt("MCINT requires two input arguments.");
	} else if (nlhs > 1) {
		mexErrMsgTxt("MCINT requires one output argument.");
	}


	/* Check the dimensions of HKLE.  HKLE can be 4 X 1 or 1 X 4. */

	m = mxGetM(HKLE_IN);
	n = mxGetN(HKLE_IN);
	if (!mxIsNumeric(HKLE_IN) || mxIsComplex(HKLE_IN) || 
		mxIsSparse(HKLE_IN)  || !mxIsDouble(HKLE_IN) ||
		(max(m,n) != 4) || (min(m,n) != 1)) {
		mexErrMsgTxt("MCINT requires that HKLE be a 4 x 1 vector.");
	}

	/* Check the dimensions of D.  D can be 12 X 1 or 1 X 12. */

	m = mxGetM(D_IN);
	n = mxGetN(D_IN);
	if (!mxIsNumeric(D_IN) || mxIsComplex(D_IN) || 
		mxIsSparse(D_IN)  || !mxIsDouble(D_IN) ||
		(max(m,n) != 12) || (min(m,n) != 1)) {
		mexErrMsgTxt("MCINT requires that D be a 12 x 1 vector.");
	}

 	/* Check the dimensions of FWHM. FWHM can be 4 X 1 or 1 X 4. */

	m = mxGetM(FWHM_IN);
	n = mxGetN(FWHM_IN);
	if (!mxIsNumeric(FWHM_IN) || mxIsComplex(FWHM_IN) || 
		mxIsSparse(FWHM_IN)  || !mxIsDouble(FWHM_IN) ||
		(max(m,n) != 4) || (min(m,n) != 1)) {
		mexErrMsgTxt("MCINT requires that FWHM be a 4 x 1 vector.");
	}

	/* Check the dimensions of B. B must be 4 X 4. */

	m = mxGetM(B_IN);
	n = mxGetN(B_IN);
	if (!mxIsNumeric(B_IN) || mxIsComplex(B_IN) || 
		mxIsSparse(B_IN)  || !mxIsDouble(B_IN) ||
		(max(m,n) != 16) || (min(m,n) != 1)) {
		mexErrMsgTxt("MCINT requires that B be a 16 x 1 vector.");
	}



 
                                          
	/* Create a matrix for the return argument */

	S_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);


	/* Assign pointers to the various parameters */

	s = mxGetPr(S_OUT);

	nmc = mxGetPr(NMC_IN);
	hkle = mxGetPr(HKLE_IN);
	d = mxGetPr(D_IN);
	fwhm = mxGetPr(FWHM_IN);
	b = mxGetPr(B_IN);


	/* Do the actual computations in a subroutine */

	mcint(s,nmc,hkle,d,fwhm,b);
	return;
}



                                          

