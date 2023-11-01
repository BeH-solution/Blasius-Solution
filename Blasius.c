/*!*****************************************************************************!
 *!         Blasius.c
 *!*****************************************************************************!
 *
 *!=============================================================================!
 *|
 *|   Program:   Blasius (double precision version)
 *|
 *|   Coding:    BeH Research Lab.
 *|   Compiler:  ANSI C / Microsoft Visual Studio
 *|
 *|   Purpose:   Solve two point boundary value problem using "shooting method"
 *|              for the "Blasius Boundary Layer Equation"
 *|
 *|              Governing Equation: (non-linear 3rd-order ODE)
 *|
 *|                    f'''(eta) + f*f''/2 = 0
 *|                    eta = y/SQRT(vx/U_inf)
 *|                    v: kinematic viscosity
 *|
 *|                (Boundary Conditions)
 *|
 *|                    f (0)   = 0.0
 *|                    f'(0)   = 0.0
 *|                    f'(inf) = 1.0
 *|
 *|                (1st order ODE sytem of Blasius equation)
 *|
 *|                    g0' =  g1           g0(0)   = 0.0
 *|                    g1' =  g2           g1(0)   = 0.0
 *|                    g2' = -g0*g2 / 2    g1(inf) = 1.0
 *|
 *|                 where g0 = f, g1 = f', g2 = f''
 *|
 *|
 *|   Supplied Functions:
 *|
 *|     function  Shoot_NR     : shooting method by Newton Rhapson
 *|     function  Shoot_Bisec  : shooting method by Bisection
 *|     function  Shoot_SysEqn : shooting method by system eqn. solver routine
 *|     function  RK4th        : 4th order Runge-Kutta method (3/8 rule)
 *|     function  RK5th        : 5th order Runge-Kutta method (Butcher's)
 *|
 *|   Notes:
 *|
 *!=============================================================================!
 */
 // #define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>                 /* for NULL & BUFSIZ   */
#include <stdlib.h>                /* for exit(1)         */
#include <malloc.h>                /* for malloc & free   */
#include <math.h>                  /* for math functions  */

#define   newline    printf("\n")


void  Shoot_NR(int, double, double, double, double*, double,
    void (*)(double, double*, double*));
void  Shoot_Bisec(int, double, double, double, double*, double,
    void (*)(double, double*, double*));
void  Shoot_SysEqn(int, double, double, double, double*, double,
    void (*)(double, double*, double*));

void  RK5th(int, void (*)(double, double*, double*), double*, double*, double);
void  RK4th(int, void (*)(double, double*, double*), double*, double*, double);
void  Write_Result(FILE*, char*, double, double*, int);

void  ludcmp(double**, int, int*);
void  lubksb(double**, int, int*, double[]);
void  pause(char[]);


/*!----------------------------------------------------------------------------!
 *|   Main Function
 *!----------------------------------------------------------------------------!
 */
int   main(int argc, char* argv[])
{
    int      ne = 3, nc = 1;
    int      i, imax, iprn, flag;
    double   a, b, eta, h, eps;
    double   f[3], f0, f1, f2;
    char* outfile, * method;


    FILE* fpout;

    void     BlasiusEqn(double, double*, double*);


    /*----------------------------------------------------------------!
     *     Set condition parameters & Open output file
     *----------------------------------------------------------------!
     */
    a = 0.00;               // low  limit of eta
    b = 12.00;              // high limit of eta
    h = 0.01;               // increment of eta value
    f[0] = 0.0;             // f(0)  = 0.0
    f[1] = 0.0;             // f'(0) = 0.0 
    f[2] = 1.0;             // f'(inf) = 1.0

    eps = 0.1e-9;
    iprn = 10;


    /*----------------------------------------------------------------!
     *     Choose method & Execute shooting to find initial values
     *----------------------------------------------------------------!
     */
    method = "Newton-Rhapson method";
    outfile = "Blasius_01.out";
    // method  = "Bisection method";
    // outfile = "Blasius_02.out";
    // method  = "System Eqn. Solver";
    // outfile = "Blasius_03.out";

    newline;
    fprintf(stderr, "      Shooting by %s\n", method);
    fprintf(stderr, "      output file : \'%s\'\n", outfile);
    newline;

    Shoot_NR(ne, a, b, h, f, eps, BlasiusEqn);
    // Shoot_Bisec(ne, a, b, h, f, eps, BlasiusEqn);
    // Shoot_SysEqn(ne, a, b, h, f, eps, BlasiusEqn);


    /*----------------------------------------------------------------!
     *     Save initial value & Write result
     *----------------------------------------------------------------!
     */
     //  fpout = fopen( outfile, "w" ); 
    fopen_s(&fpout, outfile, "w");
    eta = a;
    f0 = f[0]; f1 = f[1]; f2 = f[2];  // f(0), f'(0) and f''(0) after computing

    Write_Result(fpout, method, eta, f, flag = 1);

    imax = (int)((b - a) / h);
    for (i = 0; i < imax; i++) {

        RK5th(ne, BlasiusEqn, &eta, f, h);     //  eta = eta + h

        if ((i + 1) / iprn * iprn == i + 1)
            Write_Result(fpout, method, eta, f, flag = 0);
    }

    newline;
    fprintf(stderr, "      del. eta = %7.3f\n", h);
    fprintf(stderr, "      max. eta = %16.12f\n", eta);
    newline;
    fprintf(stderr, "      Calculated I.C's : [%4.2f, %4.2f, %16.12f]\n", f0, f1, f2);
    fprintf(stderr, "      f(inf)   = %16.12f\n", f[0]);
    fprintf(stderr, "      f'(inf)  = %16.12f\n", f[1]);
    fprintf(stderr, "      f''(inf) = %16.12f\n", f[2]);
    newline;


    fclose(fpout);
    fprintf(stderr, ">>> Execution is completed... \n");
    system("pause");
    return 0;
}



/*=============================================================================!
 *|
 *|   Function:  BlasiusEqn (double precision version)
 *|
 *|   Purpose:   Set single ODE for high order ODE system vector
 *|
 *|              Governing Equation: (non-linear 3rd-order ODE)
 *|
 *|                    f'''(eta) + f*f''/2 = 0
 *|                    eta = y/SQRT(vx/U_inf)
 *|                    v: kinematic viscosity
 *|
 *|                (Boundary Conditions)
 *|
 *|                    f (0)   = 0.0
 *|                    f'(0)   = 0.0
 *|                    f'(inf) = 1.0
 *|
 *|                (1st order ODE sytem of Blasius equation)
 *|
 *|                    f0' =  f1           f0(0)   = 0.0
 *|                    f1' =  f2           f1(0)   = 0.0
 *|                    f2' = -f0*f2 / 2    f1(inf) = 1.0
 *|
 *|                 where f0 = f, f1 = f', f2 = f''
 *|
 *|
 *|   Usage:     BlasiusEqn( x, f, dfdx );
 *|
 *|   Arguments:
 *|     x       - independent variable [i]
 *|     f[i]    - vector of dependent variable [i]
 *|     dfdx[i] - functions of Ordinary Differential Equation [o]
 *|
 *|   Reference: Cheney, E.W., and Kincaid,
 *|            - Numerical Math. and Computing
 *|
 *!============================================================================!
 */
void  BlasiusEqn(double x, double* f, double* dfdx)
{
    dfdx[0] = f[1];
    dfdx[1] = f[2];
    dfdx[2] = -f[0] * f[2] / 2.0;

    return;
}


/*=============================================================================!
 *|
 *|   Function:  Write_Result (double precision version)
 *|
 *|   Purpose:   print out the results to given file pointer
 *|
 *|   Usage:     Write_Result( fp, method, eta, f, flag );
 *|
 *|   Arguments:
 *|     fp      - file pointer [i]
 *|     method  - name of solving method [i]
 *|     eta     - independent variable [i]
 *|     f[i]    - dependent data array [i]
 *|     flag    - flag for header printing [i]
 *|
 *!============================================================================!
 */
void  Write_Result(FILE* fp, char* method, double eta, double* f, int flag)
{
    if (flag == 1) {
        fprintf(fp, "\n         The numerical solution of Blasius equation");
        fprintf(fp, "\n             shooting by %s", method);
        fprintf(fp, "\n\n\n     ");
        fprintf(fp, "eta      f(eta)        f'(eta)       f''(eta)         V    \n   ");
        fprintf(fp, "--------------------------------------------------------------\n");
    }

    fprintf(fp, "%9.4f%14.9f%14.9f%14.9f%14.9f\n",
        eta, f[0], f[1], f[2], (eta * f[1] - f[0]) / 2.0);

    return;
}


/*=============================================================================!
 *|
 *|   Function:  Shoot_NR (double precision version)
 *|
 *|   Coding:    BeH Research Lab.
 *|   Compiler:  ANSI C / Microsoft Visual Studio
 *|
 *|   Purpose:   Find unknown initial value by "shooting method"
 *|              from given boudary value with Newton-Rhapson method.
 *|
 *|   Usage:     Shoot_NR( ne,a,b,h,bv,eps,f_sys );
 *|
 *|   Arguments:
 *|     ne     - number of ODE system equations [i]
 *|     a      - starting (initial) point [i]
 *|     b      - ending  (boundary) point [i]
 *|     h      - increment of independent variable from a to b [i]
 *|     bv[i]  - array of boundary value and returned initial values [i/o]
 *|     eps    - error criterion for finding exact unknown IC [i]
 *|     f_sys  - name of subroutine which support the system equation [i]
 *|
 *|     f[i]   - array of IC vector : function of x "f(x)"
 *|     z      - unknown IC: returned as computed value
 *|              behaves as independent variable of function g(z)
 *|     g      - unknown BC function : behaves as function of z "g(z)"
 *|
 *|   Reference: Author
 *|            - reference Name
 *|
 *!============================================================================!
 */
void  Shoot_NR(int ne, double a, double b, double h, double* bv, double eps,
    void (*f_sys)(double, double*, double*))
{
    int      i, imax, iter;
    double   x, f0, f1, f2, error;
    double* f, z0, z, g0, g, dgdz;


    /*----------------------------------------------------------------!
     *     Compute 1st IC w/ guess value
     *----------------------------------------------------------------!
     */

    f = (double*)malloc((size_t)(ne * sizeof(double)));  // memory allocation

    f0 = f[0] = bv[0];               //  save f(0)
    f1 = f[1] = bv[1];               //  save f'(0)
    f2 = f[2] = 0.0;                //  f''(0) by arbitrary


    imax = (int)((b - a) / h);
    x = a;
    for (i = 0; i < imax; i++)     //  get the solution at boundary point
        RK5th(ne, f_sys, &x, f, h);

    g = error = bv[2] - f[1];       //  f[1] = f'(inf) after compute
    dgdz = 0.1;                     //  any increment at 1st time
    z0 = f2; g0 = g;                 //  substitute to old value
    z = f2 + dgdz;                   //  set unknown initial value


   /*----------------------------------------------------------------!
    *     Change initial value and check error
    *----------------------------------------------------------------!
    */
    iter = 0;
    while (fabs(error) > eps) {
        x = a;
        f[0] = f0; f[1] = f1; f[2] = z;    //  set IC with new value: z
        for (i = 0; i < imax; i++)
            RK5th(ne, f_sys, &x, f, h);
        // x = x + h

        g = error = bv[2] - f[1];         // bv[2] = f'(inf)
        dgdz = g / (g - g0) * (z - z0);
        z0 = z;  g0 = g;                   // substitute to old value
        z -= dgdz;                         // Newton's iteration
        iter++;
    }
    fprintf(stderr, "      Max. epsilon     = %8.2e\n", eps);
    fprintf(stderr, "      Newton iteration = %4d  \n", iter);


    /*----------------------------------------------------------------!
     *     Restore initial values & Clear memory allocation
     *----------------------------------------------------------------!
     */
    bv[0] = f0;  bv[1] = f1; bv[2] = z;

    free((char*)f);
    return;
}


/*!============================================================================!
 *|
 *|   Function:  Shoot_Bisec (double precision version)
 *|
 *|   Coding:    BeH Research Lab.
 *|   Compiler:  ANSI C / Microsoft Visual Studio
 *|
 *|   Purpose:   Find unknown initial value by "shooting method"
 *|              from given boudary value with bisection method.
 *|
 *|   Usage:     Shoot_Bisec( ne,a,b,h,bv,eps,f_sys );
 *|
 *|   Arguments:
 *|     ne     - number of ODE system equations
 *|     a      - starting (initial) point [i]
 *|     b      - ending  (boundary) point [i]
 *|     h      - increment of independent variable from a to b [i]
 *|     bv[i]  - array of boundary value and returned initial values [i/o]
 *|     eps    - error criterion for finding exact unknown IC [i]
 *|     f_sys  - name of subroutine which support the system equation [i]
 *|
 *|     nc     - number of unknown IC (starting values)
 *|     f(i)   - array of IC vector [i/o] : function of x "f(x)"
 *|     z(i)   - array of unknown IC vector
 *|              behaves as independent variable of function g(z)
 *|     g(i)   - array of   known BC vector [i]
 *|              behaves as function of z "g(z)"
 *|
 *|   Reference: Author
 *|            - Reference Name
 *|
 *!============================================================================!
 */
void  Shoot_Bisec(int ne, double a, double b, double h, double* bv, double eps,
    void (*f_sys)(double, double*, double*))
{
    int      i, imax, iter, maxit = 200;
    double* f, f0, f1, f2;
    double   x, z1, z2, g1, g2, zz, dz, zm;

    /*----------------------------------------------------------------!
     *     Compute 2 initial point values of function g(z)
     *----------------------------------------------------------------!
     */

    f = (double*)malloc((size_t)(ne * sizeof(double)));  // memory allocation

    f0 = f[0] = bv[0];               //  save f(0)
    f1 = f[1] = bv[1];               //  save f'(0)
    f2 = f[2] = 0.0;                 //  f''(0) by arbitrary 1

    imax = (int)((b - a) / h);
    x = a;
    z1 = f[2];
    for (i = 0; i < imax; i++)     //  get the solution at boundary point
        RK5th(ne, f_sys, &x, f, h);

    g1 = bv[2] - f[1];               //  f[1] = f'(inf) after compute


    f2 = f[2] = 1.0;                 //  f''(0) by arbitrary 2
    x = a;
    z2 = f[2];
    for (i = 0; i < imax; i++)     //  get the solution at boundary point
        RK5th(ne, f_sys, &x, f, h);

    g2 = bv[2] - f[1];               //  f[1] = f'(inf) after compute

    if (g1 * g2 >= 0.0) pause("*** No Bracket, Set IC again...");

    /*----------------------------------------------------------------!
     *     Start bisection
     *----------------------------------------------------------------!
     */
    if (g1 < 0.0) {
        zz = z1;  dz = z2 - z1;
    }
    else {
        zz = z2;  dz = z1 - z2;
    }

    iter = 0;
    while (fabs(dz) > eps / 100.0) {
        iter++;
        if (iter > maxit) {
            pause("*** Too many bisections in Root_Bisec...");
            exit(0);
        }
        dz = dz / 2.0;
        zm = zz + dz;

        x = a;
        f[0] = f0;
        f[1] = f1;
        f[2] = zm;                 // new bisectered value
        for (i = 1; i < imax; i++)
            RK5th(ne, f_sys, &x, f, h);

        g2 = bv[2] - f[1];
        if (g2 <= 0.0) zz = zm;
    }
    fprintf(stderr, "      Max. epsilon       = %8.2e\n", eps);
    fprintf(stderr, "      Bisction iteration = %4d  \n", iter);


    /*----------------------------------------------------------------!
     *     Restore initial values & Clear memory allocation
     *----------------------------------------------------------------!
     */
    bv[0] = f0;  bv[1] = f1; bv[2] = zz;

    free((char*)f);
    return;
}


/*=============================================================================!
 *|
 *|   Function:  Shoot_SysEqn (double precision version)
 *|
 *|   Coding:    BeH Research Lab.
 *|   Compiler:  ANSI C / Microsoft Visual Studio
 *|
 *|   Purpose:   Find unknown initial value by "shooting method"
 *|              from given boudary value with system eqn, solver routine.
 *|
 *|   Usage:     Shoot_SysEqn( ne,a,b,h,bv,eps,f_sys );
 *|
 *|   Arguments:
 *|     ne     - number of ODE system equations [i]
 *|     a      - starting (initial) point [i]
 *|     b      - ending  (boundary) point [i]
 *|     h      - increment of independent variable from a to b [i]
 *|     bv[i]  - array of boundary value and returned initial values [i/o]
 *|     eps    - error criterion for finding exact unknown IC [i]
 *|     f_sys  - name of subroutine which support the system equation [i]
 *|
 *|     nc     - number of unknown IC (starting values)
 *|     f(i)   - array of IC vector
 *|     z(i)   - array of unknown IC vector
 *|     delz   - increment of unknown starting value vector
 *|     dz(nc) - value of discrepancy F at adjustable parameter z+delz
 *|     ajparm - adjust parameter (set exactly through trial & error)
 *|
 *|   Reference: W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling
 *|            - "Numerical Recipes", Cambridge Univ. Press. 1986
 *|
 *!============================================================================!
 */
void  Shoot_SysEqn(int ne, double a, double b, double h, double* bv, double eps,
    void (*f_sys)(double, double*, double*))
{
    int      i, j, imax, indx, nc = 1;
    double   sav, x, f0, f1, f2, * f;
    double* z, * delz, * dz, * dz_old, ** dfdz, ajparm;


    /*----------------------------------------------------------------!
     *     Memory allocation for array
     *----------------------------------------------------------------!
     */
    f = (double*)malloc((size_t)(ne * sizeof(double)));

    z = (double*)malloc((size_t)(nc * sizeof(double)));
    delz = (double*)malloc((size_t)(nc * sizeof(double)));
    dz = (double*)malloc((size_t)(nc * sizeof(double)));
    dz_old = (double*)malloc((size_t)(nc * sizeof(double)));

    dfdz = (double**)malloc((size_t)(nc * sizeof(double*)));
    for (i = 0; i < nc; i++)
        dfdz[i] = (double*)malloc((size_t)(nc * sizeof(double)));



    /*----------------------------------------------------------------!
     *     Compute new boundary value
     *----------------------------------------------------------------!
     */
    imax = (int)((b - a) / h);

    z[0] = 0.0;                      //  f''(0) by arbitrary
    delz[0] = 0.01;                  //  increment of unknown I.C's

    // set I.C's
    f0 = f[0] = bv[0];               //  save f(0)
    f1 = f[1] = bv[1];               //  save f'(0)
    f2 = f[2] = z[0];                //  save f''(0)

    x = a;
    for (i = 0; i < imax; i++)  RK5th(ne, f_sys, &x, f, h);

    // Set B.C's  : dz_old is value of discrepancy F at adjustable parameter z+delz
    ajparm = 3.00717734818;          // ajparm = 3.21224629870;
    dz_old[0] = f[1] - ajparm;

    for (j = 0; j < nc; j++) {
        sav = z[j];
        z[j] += delz[j];           // increment parameter j

        // set I.C's
        f[0] = bv[0];              //  save f(0)
        f[1] = bv[1];              //  save f'(0)
        f[2] = z[j];               //  save f''(0) again

        x = a;
        for (i = 0; i < imax; i++)  RK5th(ne, f_sys, &x, f, h);

        dz[j] = f[1] - ajparm;     // set boundary values
        for (i = 0; i < nc; i++)  dfdz[i][j] = (dz[i] - dz_old[i]) / delz[j];

        z[j] = sav;
    }
    for (j = 0; j < nc; j++)  dz[j] = -dz_old[j];

    //-----------------------------------!
    ludcmp(dfdz, nc, &indx);
    lubksb(dfdz, nc, &indx, dz);
    //-----------------------------------!

    for (j = 0; j < nc; j++)   z[j] += dz[j];

    fprintf(stderr, "      Max. epsilon = %8.2e\n", eps);

    /*----------------------------------------------------------------!
     *     Restore initial values & Clear memory allocation
     *----------------------------------------------------------------!
     */
    bv[0] = f0;  bv[1] = f1; bv[2] = z[0];

    for (i = 0; i < nc; i++)  free((char*)dfdz[i]);
    free((char*)dfdz);
    free((char*)dz_old);
    free((char*)dz);
    free((char*)delz);
    free((char*)z);
    free((char*)f);

    return;
}


/*=============================================================================!
 *|
 *|   Function:  RK5th (double precision version)
 *|
 *|   Purpose:   5th order Runge-Kutta-Butcher's method to solve ODE system
 *|
 *|   Usage:     RK5th( ne,f_sys,&x,y,h );
 *|
 *|   Arguments:
 *|     ne     - number of ODE system equation [i]
 *|     f_sys  - subroutine for ODE system     [i]
 *|     x      - ODE independent variable      [i/o]
 *|              it is replaced to next step
 *|     y      - ODE dependent variable vector [i/o]
 *|     h      - step size of increment        [i]
 *|  // pole   - status for the solution approach to a stiff pole [o]
 *|
 *|   Required:  Call system equation subroutine "f_sys" in order to input
 *|              ODE system
 *|
 *|   Reference: 1. Curtis F. Gerald, Patrick O. Wheatley,
 *|                "Applied Numerical Analysis", 4e, Addison Wesley, (1989)
 *|              2. S.C.Chapra, R.P.Canale,
 *|                "Numerical Methods for Engineers", 4e, (2003)
 *|
 *!============================================================================!
 */
void  RK5th(int ne, void (*f_sys)(double, double*, double*),
    double* x, double* y, double h)
{
    int      k;
    double   hh, xparm, * yparm;
    double* RK1, * RK2, * RK3, * RK4, * RK5, * RK6;


    /*----------------------------------------------------------------!
     *     Memory allocation for array
     *----------------------------------------------------------------!
     */
    yparm = (double*)malloc((size_t)(ne * sizeof(double)));
    RK1 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK2 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK3 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK4 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK5 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK6 = (double*)malloc((size_t)(ne * sizeof(double)));


    /*----------------------------------------------------------------!
     *     Start procedure
     *----------------------------------------------------------------!
     */

    (*f_sys)(*x, y, RK1);        //----------------------!  Compute k1

    hh = h / 4.0;                  //----------------------!  Compute k2
    xparm = *x + hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * RK1[k];

    (*f_sys)(xparm, yparm, RK2);

    hh = h / 8.0;                  //----------------------!  Compute k3
    xparm = *x + 2.0 * hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * (RK1[k] + RK2[k]);

    (*f_sys)(xparm, yparm, RK3);

    hh = h / 2.0;                  //----------------------!  Compute k4
    xparm = *x + hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * (-RK2[k] + 2.0 * RK3[k]);

    (*f_sys)(xparm, yparm, RK4);

    hh = h / 16.0;                 //----------------------!  Compute k5
    xparm = *x + 8.0 * hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * (3.0 * RK1[k] + 9.0 * RK4[k]);
    //  yparm[k] = y[k] + hh*( 3.0*RK1[k] +10.0*RK4[k] );

    (*f_sys)(xparm, yparm, RK5);

    hh = h / 7.0;                  //----------------------!  Compute k5
    xparm = *x + h;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * (-3.0 * RK1[k] + 2.0 * RK2[k]
            + 12.0 * RK3[k] - 12.0 * RK4[k]
            + 8.0 * RK5[k]);

    (*f_sys)(xparm, yparm, RK6);

    hh = h / 90.0;                 //----------------------!  sum it
    *x += h;
    for (k = 0; k < ne; k++)
        y[k] = y[k] + hh * (7.0 * RK1[k] + 32.0 * RK3[k]
            + 7.0 * RK6[k] + 32.0 * RK5[k]
            + 12.0 * RK4[k]);


    /*----------------------------------------------------------------!
     *     Clear memory allocation
     *----------------------------------------------------------------!
     */
    free((char*)RK6);
    free((char*)RK5);
    free((char*)RK4);
    free((char*)RK3);
    free((char*)RK2);
    free((char*)RK1);
    free((char*)yparm);

    return;
}


/*!============================================================================!
 *|
 *|   Purpose:   3/8 4th order Runge-Kutta method to solve ODE system
 *|
 *|   Usage:     RK4th( ne,f_sys,&x,y,h );
 *|
 *|   Arguments:
 *|     ne     - number of ODE system equation [i]
 *|     f_sys  - subroutine for ODE system     [i]
 *|     x      - ODE independent variable      [i/o]
 *|              it is replaced to next step
 *|     y      - ODE dependent variable vector [i/o]
 *|     h      - step size of increment        [i]
 *|  // pole   - status for the solution approach to a stiff pole [o]
 *|
 *|   Required:  Call system equation subroutine "f_sys" in order to input
 *|              ODE system
 *|
 *|   Reference: 1. B.Carnahan, H.A.Luther, J.O.Wilkes,
 *|                "Applied Numerical Methods", Wiley (1969)
 *|
 *!============================================================================!
 */
void  RK4th(int ne, void (*f_sys)(double, double*, double*),
    double* x, double* y, double h)
{
    int      k;
    double   hh, xparm, * yparm;
    double* RK1, * RK2, * RK3, * RK4;


    /*----------------------------------------------------------------!
     *     Memory allocation for array
     *----------------------------------------------------------------!
     */
    yparm = (double*)malloc((size_t)(ne * sizeof(double)));
    RK1 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK2 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK3 = (double*)malloc((size_t)(ne * sizeof(double)));
    RK4 = (double*)malloc((size_t)(ne * sizeof(double)));


    /*----------------------------------------------------------------!
     *     Start procedure
     *----------------------------------------------------------------!
     */

    (*f_sys)(*x, y, RK1);        //----------------------!  Compute k1

    hh = h / 3.0;                  //----------------------!  Compute k2
    xparm = *x + hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * RK1[k];

    (*f_sys)(xparm, yparm, RK2);

    hh = h / 3.0;                  //----------------------!  Compute k3
    xparm = *x + 2.0 * hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * (-RK1[k] + 3.0 * RK2[k]);

    (*f_sys)(xparm, yparm, RK3);

    hh = h;                      //----------------------!  Compute k4
    xparm = *x + hh;
    for (k = 0; k < ne; k++)
        yparm[k] = y[k] + hh * (RK1[k] - RK2[k] + RK3[k]);

    (*f_sys)(xparm, yparm, RK4);

    hh = h / 8.0;                  //----------------------!  sum it
    *x += h;
    for (k = 0; k < ne; k++)
        y[k] = y[k] + hh * (RK1[k] + 3.0 * RK2[k] + 3.0 * RK3[k] + RK4[k]);


    /*----------------------------------------------------------------!
     *     Clear memory allocation
     *----------------------------------------------------------------!
     */
    free((char*)RK4);
    free((char*)RK3);
    free((char*)RK2);
    free((char*)RK1);
    free((char*)yparm);

    return;
}



/*!============================================================================!
 *|
 *|   Function: ludcmp (double precision version)
 *|
 *|   Purpose:  Decompose the matrix A as LU of a rowwise permuation
 *|             by partial pivoting.
 *|
 *|   Usage:    ludcmp( A,n,indx );
 *|
 *|   Arguments:
 *|     A      - matrix [i]
 *|     n      - n*n matrix size [i]
 *|     indx   - output vector which records the row permuation
 *|              effected by the partial pivoting [o]
 *|     d      - +1 if the number of row interchanges are even
 *|              -1 if the number of row interchanges are odd [o]
 *|
 *|   Reference: 1. W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling,
 *|                "Numerical Recipes", 2nd ed., Cambridge (1992).
 *|
 *!============================================================================!
 */
void  ludcmp(double** A, int n, int* indx)
{
    int     i, j, k, imax;
    double  d, tiny = 1.0e-20;
    double  Amax, dum, sum, * vv;

    vv = (double*)malloc((size_t)(n * sizeof(double)));

    d = 1.0;
    for (i = 0; i < n; i++) {
        Amax = 0.0;
        for (j = 0; j < n; j++)
            if (fabs(A[i][j]) > Amax)  Amax = fabs(A[i][j]);
        if (Amax == 0.0) pause("singular matrix in ludcmp");
        vv[i] = 1.0 / Amax;
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < j - 1; j++) {
            sum = A[i][j];
            for (k = 0; k < i - 1; k++)  sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
        }

        Amax = 0.0;
        for (i = j; i < n; i++) {
            sum = A[i][j];
            for (k = 0; k < j - 1; k++)  sum -= A[i][k] * A[k][j];
            A[i][j] = sum;

            dum = vv[i] * fabs(sum);
            if (dum >= Amax) {
                imax = i;  Amax = dum;
            }
        }

        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = A[imax][k];
                A[imax][k] = A[j][k];
                A[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if (A[j][j] == 0.0) A[j][j] = tiny;
        if (j != n - 1) {
            dum = 1.0 / A[j][j];
            for (i = j + 1; i < n; i++) A[i][j] *= dum;
        }
    }

    free((char*)vv);
    return;
}



/*!============================================================================!
 *|
 *|   Function: lubksb (double precision version)
 *|
 *|   Purpose:  Decompose the matrix A as LU of a rowwise permuation
 *|             by partial pivoting.
 *|
 *|   Usage:    solves the set of N linear equations A*x=B.
 *|
 *|   Arguments:
 *|     A      - matrix as LU decomposition determined by routine LUDCMP [i]
 *|     n      - n*n matrix size [i]
 *|     np     - max. physical dimension of matrix A.
 *|     indx   - the permuation vector returned by routine LUDCMP [i]
 *|     B      - column vector as input and returned with solution vector x
 *|
 *|   Reference: 1. W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling,
 *|                "Numerical Recipes", 2nd ed., Cambridge (1992).
 *|
 *!============================================================================!
 */
void  lubksb(double** A, int n, int* indx, double B[])
{
    int      i, j, ii = 0, ll;
    double   sum;

    for (i = 0; i < n; i++) {
        ll = indx[i];
        sum = B[ll];
        B[ll] = B[i];
        if (ii != 0)
            for (j = ii; j <= i - 1; j++) sum -= A[i][j] * B[j];
        else if (sum != 0.0)  ii = i;

        B[i] = sum;
    }

    for (i = n - 1; i >= 0; i--) {
        sum = B[i];
        for (j = i + 1; j < n; j++) sum -= A[i][j] * B[j];
        B[i] = sum / A[i][i];
    }

    return;
}



/*!----------------------------------------------------------------------------!
 *|         Implementation of pause function
 *!----------------------------------------------------------------------------!
 */
void  pause(char message[])
{
    char ch;

    if (message == "")  message = "Pause: Press any key to continue...";
    fprintf(stderr, "\n%s", message);
    while (1)                    /* endless loop */
        if (ch = getchar() == '\n') break;

    return;
}
