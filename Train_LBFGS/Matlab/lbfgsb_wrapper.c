/*=================================================================
 *
 * lbfgsb_wrapper.c, mex gatewate interface
 * You can call this directly from Matlab,
 * but I strongly suggest using lbfgsb.m, since many
 * variables are modified in-place, and lbfgsb.m will
 * handle them correctly. It also shows you how to use
 * this routine
 *  (so just modify lbfgsb.m to suit your own taste)
 *
 * Written by Stephen Becker, Feb 14 2012
 * Major update Feb 21 2015
 *
 * This re-distributes a C version of lbfgsb 3.0
 * (free; see http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html)
 *
 * There are also versions 2.1 and 2.4 of the lbfgsb library
 *      For v 2.1, Peter Carbonetto's mex interface works; see
 * http://www.mathworks.com/matlabcentral/fileexchange/15061-matlab-interface-for-l-bfgs-b
 *      and also
 *  http://www.cs.ubc.ca/~pcarbo/lbfgsb-for-matlab.html
 *
 *      For v 2.4 (I don't know where you can find this version though),
 *      use the mex files from here:
 * http://www.cs.toronto.edu/~liam/software.shtml
 *
 *  
 *  You should use version 3.0! It's better.  The old versions are basically
 *  15 years old; changes in version 3.0 are fundamental (see the 2011 "remark" paper
 *  listed at Nocedal's website)
 *
 *  Also, this version has converted everything to C, so there is no mixed
 *  C and fortran, and compilation is much easier. In the conversion,
 *  I changed a few conventions, and some of these are marked with my 
 *  initials SRB
 *      -- Stephen Becker, Feb 2015. stephen.becker@colorado.edu
 *
 *
 *
 * Inputs (in order):
 *  m       # of memory vectors to keep
 *  x       initial guess or current point, also used to determine the size of the problem (Nx1)
 *  l       list of upper bounds (Nx1)
 *  u       list of lower bounds (Nx1)
 *  nbd     list of which bounds are active (Nx1):
 *              0=neither u nor l, 1 = only l, 2 = both u and l, 3 = only u
 *  f       value of function (1x1)
 *  g       gradient of function (Nx1)
 *
 *  factr   stopping crit: 1e+12 for low accuracy, 1e7 for moderate, 1e1 for high accuracy
 *              (will be multiplied by machine precision)
 *  pgtol   stopping crit for infinity norm of projected gradient
 *  wa      work space array (double)
 *  iwa     work space array (int)
 *  task    string                                  SRB UPDATE: now an integer
 *  iprint  how verbose the fortran program should be
 *  csave       some iformation (string, length 60) SRB UPDATE: now an integer
 *  lsave[4]    some information (logicals)
 *  isave[44]   some information (integers)
 *  dsave[29]   some information (doubles)
 *
 * Outputs:
 *  f, task, csave, lsave, isave, dsave
 *
 * Warning: the following variables are modified in-place
 *  x, g, wa, iwa
 *
 *=================================================================*/
#include <math.h>
/* #include "mex.h" */  /* now, mex.h included in lbfgsb.h */

#include "lbfgsb.h"

#include <string.h>
#include <limits.h> /* for CHAR_BIT */
#include <assert.h>


/* These are the fortran arguments in order */
#define N_m 0
#define N_x 1
#define N_l 2
#define N_u 3
#define N_nbd 4
#define N_factr 6
#define N_pgtol 7
#define N_iprint 8
#define N_iterMax 9 /* new */
#define N_total_iterMax 10 /* new */

#define N_fcn 5 /* new in 2015 */

/* these depend on the lbfgsb release version */
#define LENGTH_STRING 60
#define LENGTH_LSAVE 4
#define LENGTH_ISAVE 44
#define LENGTH_DSAVE 29


#ifdef DEBUG
#define debugPrintf mexPrintf
#else
#define debugPrintf fakePrintf
#endif
/* If the 'DEBUG' symbol is undefined, then don't print: */
int fakePrintf(const char *format, ...){
    return 0;
}

mxLogical isInt( const mxArray *pm ) {
    /* We typedef "integer" to be "long int", and this is not
     * constant across different computers.
     *
     * CHAR_BIT is from limits.h
     * If using gcc, you can run `gcc -dM -E - < /dev/null | grep CHAR_BIT`
     *  and it should define the symbol __CHAR_BIT__, so this is another way.
     *
     *  This will match the typdef "integer" which is what
     *  the lbfgsb codes uses for integers.
     * */
    
    /* debugPrintf("Sizeof(int) is %d\n", sizeof(int) ); */
    switch ( CHAR_BIT * sizeof(integer) ) {
        case 16 :
            return mxIsInt16(pm);
        case 32:
            return mxIsInt32(pm);
        case 64:
            return mxIsInt64(pm);
        default:
            mexErrMsgTxt("You have a weird computer that I don't know how to support");
            return false;
    }
}


/* Main mex gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )   { 
    
    integer iprint = (integer)1;
    integer task=(integer)START, csave=(integer)1;
    integer iterations = 0;
    integer total_iterations = 0;

    int  iterMax = 100;
    int  total_iterMax = 200;


    integer   n, m, *nbd=NULL, *iwa=NULL; 
    double  f=0, factr, pgtol, *x, *l, *u, *g, *wa=NULL;
    int     i;
    mxLogical FREE_nbd=false;

    int ndim = 2; /* for lcc compiler, must declare these here, not later ... */
    mwSize dims[2] = { LENGTH_ISAVE, 1 };

    logical lsave[LENGTH_LSAVE];
    integer isave[LENGTH_ISAVE];
    double  dsave[LENGTH_DSAVE];
    
    double *nbd_dbl=NULL;
    long long *nbd_long=NULL;
    
    mxArray *LHS[2];
    mxArray *RHS[3];
    double *tempX, *tempG, *tempIter;
    
    /* Parse inputs. Quite boring */
    
    if (nrhs < 5 ) mexErrMsgTxt("Needs at least 5 input arguments");
    m       = (int)*mxGetPr( prhs[N_m] );
    n       = (integer)mxGetM( prhs[N_x] );
    if ( mxGetN(prhs[N_x]) != 1 ) mexErrMsgTxt("x must be a column vector");
    if ( mxGetM(prhs[N_l]) != n ) mexErrMsgTxt("l must have same size as x");
    if ( mxGetM(prhs[N_u]) != n ) mexErrMsgTxt("u must have same size as x");
    if ( mxGetM(prhs[N_nbd]) != n ) mexErrMsgTxt("nbd must have same size as x");


    if (nlhs < 2 )  mexErrMsgTxt("Should have 2 or 3 output arguments");
    if (!mxIsDouble(prhs[N_x]))
            mexErrMsgTxt("x should be of type double!\n");
    plhs[1] = mxDuplicateArray( prhs[N_x] );
    x       = mxGetPr( plhs[1] );


    l       = mxGetPr( prhs[N_l] );
    u       = mxGetPr( prhs[N_u] );
    if ( isInt( prhs[N_nbd] ) ) {
        nbd     = (integer *)mxGetData( prhs[N_nbd] ); 
    } else {
        debugPrintf("Converting nbd array to integers\n" );
        if (!mxIsDouble(prhs[N_nbd])){
            if (mxIsInt64(prhs[N_nbd])){
                nbd_long = mxGetData( prhs[N_nbd] );
                nbd     = (integer *)mxMalloc( n * sizeof(integer) );
                assert( nbd != NULL );
                FREE_nbd = true;
                /* convert nbd_dbl (in double format) to integers */
                for (i=0;i<n;i++)
                    nbd[i]  = (integer)nbd_long[i];
            } else {
                debugPrintf("Sizeof(int) is %d bits, sizeof(integer) is %d bits\n",
                        CHAR_BIT*sizeof(int),CHAR_BIT*sizeof(integer) );
                /* integer is aliased to 'long int' and should be at least
                 * 32 bits. 'long long' should be at least 64 bits.
                 * On 64-bit Windows, it seems 'long int' is exactly 32 bits,
                 * while on 64-bit linux and Mac, it is 67 bits */
                debugPrintf("Nbd is of type %s\n", mxGetClassName( prhs[N_nbd] ) );
                mexErrMsgTxt("Nbd array not doubles or type int64!\n");
            }
        } else {
            nbd_dbl = mxGetPr( prhs[N_nbd] );
            nbd     = (integer *)mxMalloc( n * sizeof(integer) );
            assert( nbd != NULL );
            FREE_nbd = true;
            /* convert nbd_dbl (in double format) to integers */
            for (i=0;i<n;i++)
                nbd[i]  = (integer)nbd_dbl[i];
        }
    }


    /* some scalar parameters */
    if ( nrhs < N_factr+1 ) 
        factr   = 1.0e7;
    else if (mxGetNumberOfElements( prhs[N_factr] )!=1)
        factr   = 1.0e7;
    else {
        factr   = (double)mxGetScalar( prhs[N_factr] );
        if (factr < 0 )
            mexErrMsgTxt("factr must be >= 0\n");
    }

    if ( nrhs < N_pgtol+1 ) 
        pgtol   = 1.0e-5;
    else if (mxGetNumberOfElements( prhs[N_pgtol] )!=1)
        pgtol   = 1.0e-5;
    else {
        pgtol   = (double)mxGetScalar( prhs[N_pgtol] );
        if (pgtol < 0)
            mexErrMsgTxt("pgtol must be >= 0\n");
    }
    if ( nrhs < N_iprint+1 ) {
        iprint  = (integer)1;
    } else if (mxGetNumberOfElements( prhs[N_iprint] )!=1) {
        iprint  = (integer)1;
    } else {
        iprint = (integer)mxGetScalar( prhs[N_iprint] );
    }
    
    if ( nrhs >= N_iterMax+1 ) 
        iterMax = (int)mxGetScalar( prhs[N_iterMax] );
    if ( nrhs >= N_total_iterMax+1 ) 
        total_iterMax = (int)mxGetScalar( prhs[N_total_iterMax] );
    
    /* allocate memory for arrays */
    g   = (double *)mxMalloc( n * sizeof(double) );
    assert( g != NULL );
    wa      = (double *)mxMalloc( (2*m*n + 5*n + 11*m*m + 8*m ) * sizeof(double) );
    assert( wa != NULL );
    iwa     = (integer *)mxMalloc( (3*n)*sizeof(integer) );
    assert( iwa != NULL );
    

            
    /* -- Finally, done with parsing inputs. Now, call lbfgsb fortran routine */
    
    /* Be careful! This modifies many variables in-place! 
     * Basically, anything without a '&' before it will be changed in the Matlab
     * workspace */
    
    if ( nrhs < N_fcn - 1 )
        mexErrMsgTxt("For this f(x) feature, need more input aguments\n");
    RHS[0] = mxDuplicateArray( prhs[N_fcn] );
    RHS[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    RHS[2] = mxCreateDoubleScalar( 0.0 ); /* The iterations counter */
    tempX = (double*)mxGetPr( RHS[1] );
    if (!mxIsDouble(RHS[2]))
        mexErrMsgTxt("Error trying to create RHS[2]\n");
    tempIter = (double*)mxGetPr( RHS[2] );

    while ( (iterations < iterMax) && (total_iterations < total_iterMax ) ){
        total_iterations++;

        setulb(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,&task,&iprint,
                &csave,lsave,isave,dsave); /* (ftnlen) TASK_LEN, (ftnlen) CSAVE_LEN); */


        if ( IS_FG(task) ) {

            /* copy data from x to RHS[1] or just set pointer with mxSetPr */
            for (i=0;i<n;i++)
                tempX[i] = x[i];
            /*Try being bold: */
            /*mxSetPr( RHS[1], x ); */

            *tempIter = (double)iterations;
            mexCallMATLAB(2,LHS,3,RHS,"feval");
            f = mxGetScalar( LHS[0] );
            if (mxGetM(LHS[1]) != n )
                mexErrMsgTxt("Error with [f,g]=fcn(x) : g wrong size\n");
            if (mxGetN(LHS[1]) != 1 )
                mexErrMsgTxt("Error with [f,g]=fcn(x) : g wrong size (should be column vector)\n");

            /* could use memcpy, or just do it by hand... */
            if (!mxIsDouble(LHS[1]))
                mexErrMsgTxt("[f,g]=fcn(x) did not return g as type double\n");
            tempG = mxGetPr( LHS[1] );
            for (i=0;i<n;i++)
                g[i] = tempG[i];
            /* Or, be a bit bolder: */
            /*g = tempG; // Hmm, crashed */

            continue;
        }
        if ( task==NEW_X ) {
            iterations++;
            continue;
        } else
            break;

    }

    mxDestroyArray( LHS[0] );
    mxDestroyArray( LHS[1] );
    mxDestroyArray( RHS[0] );
    mxDestroyArray( RHS[1] );
            

    
    plhs[0] = mxCreateDoubleScalar( f );
    if ( nlhs >= 3 )
        plhs[2] = mxCreateDoubleScalar( task );
    if ( nlhs >= 4 )
        plhs[3] = mxCreateDoubleScalar( iterations );
    if ( nlhs >= 5 )
        plhs[4] = mxCreateDoubleScalar( total_iterations );
    if ( nlhs >= 6 )
        mexErrMsgTxt("Did not expect more than 5 outputs\n");

    if (FREE_nbd)
        mxFree(nbd);
    mxFree(g);
    mxFree(wa);
    mxFree(iwa);

    return;
}
