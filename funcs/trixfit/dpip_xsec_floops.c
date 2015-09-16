/*===========================================================
 * dpip_xsec_floops.c : C function to perform the for-loops in the
 * calculation of the cross-section
 * The calling syntax is:
 *
 *	[ik,ie] = dpip_xsec_floops(Lk,Lee,Qhh,ee)
#include <math.h>
#include "mex.h"

/* Input arguments */

/*#define LK	prhs[0]*/
/*#define LEE	prhs[1]*/
/*#define LQHH	prhs[2]*/
/*#define LW	prhs[3]*/
/*#define KK	prhs[4]*/
/*#define QHH	prhs[5]*/
/*#define EE	prhs[6]*/
/*#define W	prhs[7]*/

/* Output arguments */

/*#define IK	plhs[0]*/
/*#define IE	plhs[1]*/


/* the computational routine of the mex-file */
/*static void dpip_xsec_floops(void)*/
/*	{*/

/*	}*/
/* the gateway routine */
void mexFunction(int nrhs, const mxArray*prhs[])
	{
/*	int 	*Lk	= mxGetPr(LK);*/ /* length of k */
/*	int 	*Lee	= mxGetPr(LEE);*/ /* length of ee */
/*	int	*LQhh	= mxGetPr(LQHH);*/ /* length of Qhh */
/*	int	*Lw	= mxGetPr(LW);*/ /* length of w */
/*	double 	*kk	= mxGetPr(KK);*/ 
/*	double 	*Qhh	= mxGetPr(QHH);*/
/*	double 	*w	= mxGetPr(W);	*/
/*	double 	*EE	= mxGetPr(EE);*/
	
/*	double[*LQhh] dnew;*/
/*	double[*LQhh] ik;*/
/*	double[*LQhh] dik;*/

	for (int j=1; j < 10 ; j++) {
		for (int k=1 ; k < 10)
			printf("Hallo");
		/*	dnew[k]=abs(*w[k]-*kk[j]);*/
	}
		



	}


