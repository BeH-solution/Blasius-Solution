#!*****************************************************************************!
#!         Blasius.py
#!*****************************************************************************!
#
#!=============================================================================!
#|
#|   Program:   Blasius (double precision version)
#|
#|   Coding:    BeH Research Lab.
#|   Compiler:  Python 3.0 over
#|
#|   Purpose:   Solve two point boundary value problem using "shooting method"
#|              for the "Blasius Boundary Layer Equation"
#|
#|              Governing Equation: (non-linear 3rd-order ODE)
#|
#|                    f'''(eta) + f*f''/2 = 0
#|                    eta = y/SQRT(vx/U_inf)
#|                    v: kinematic viscosity
#|
#|                (Boundary Conditions)
#|
#|                    f (0)   = 0.0
#|                    f'(0)   = 0.0
#|                    f'(inf) = 1.0
#|
#|                (1st order ODE sytem of Blasius equation)
#|
#|                    g0' =  g1           g0(0)   = 0.0
#|                    g1' =  g2           g1(0)   = 0.0
#|                    g2' = -g0*g2 / 2    g1(inf) = 1.0
#|
#|                 where g0 = f, g1 = f', g2 = f''
#|
#|
#|   Supplied Functions:
#|
#|     function  Shoot_NR     : shooting method by Newton Rhapson
#|     function  Shoot_Bisec  : shooting method by Bisection
#|     function  RK4th        : 4th order Runge-Kutta method (3/8 rule)
#|     function  RK5th        : 5th order Runge-Kutta method (Butcher's)
#|
#|   Notes:
#|
#!=============================================================================!
###


#!=============================================================================!
#|
#|   Function:  BlasiusEqn
#|
#|   Purpose:   Set single ODE for high order ODE system vector
#|
#|              Governing Equation: (non-linear 3rd-order ODE)
#|
#|                    f'''(eta) + f*f''/2 = 0
#|                    eta = y/SQRT(vx/U_inf)
#|                    v: kinematic viscosity
#|
#|                (Boundary Conditions)
#|
#|                    f (0)   = 0.0
#|                    f'(0)   = 0.0
#|                    f'(inf) = 1.0
#|
#|                (1st order ODE sytem of Blasius equation)
#|
#|                    f0' =  f1           f0(0)   = 0.0
#|                    f1' =  f2           f1(0)   = 0.0
#|                    f2' = -f0*f2 / 2    f1(inf) = 1.0
#|
#|                 where f0 = f, f1 = f', f2 = f''
#|
#|
#|   Usage:     BlasiusEqn( x, f )
#|
#|   Arguments:
#|     x       - independent variable [i]
#|     f[i]    - vector of dependent variable [i]
#|  
#|   Return value:
#|     dfdx[i] - functions of ODE [o]
#|
#|   Reference: Cheney, E.W., and Kincaid,
#|            - Numerical Math. and Computing
#|
#!===========================================================================!
###
def BlasiusEqn( x, f ):

    dfdx = [ f[1], f[2],  -f[0]*f[2]/2.0 ]
    
    return dfdx


#!============================================================================!
#|
#|   Purpose:   3/8 4th order Runge-Kutta method to solve ODE system
#|
#|   Usage:     RK4th( ne, f_sys, x, y, h )
#|
#|   Arguments:
#|     ne     - number of ODE system equation [i]
#|     f_sys  - subroutine for ODE system     [i]
#|     x      - ODE independent variable      [i/o]
#|              it is replaced to next step
#|     y      - ODE dependent variable vector [i/o]
#|     h      - step size of increment        [i]
#|  // pole   - status for the solution approach to a stiff pole [o]
#|
#|   Required:  Call system equation subroutine "f_sys" in order to input
#|              ODE system
#|
#|   Reference: 1. B.Carnahan, H.A.Luther, J.O.Wilkes,
#|                "Applied Numerical Methods", Wiley (1969)
#|
#!============================================================================!
#
def RK4th( ne, f_sys, x, y, h ):

    #----------------------------------------------------------------!
    #     Start procedure
    #----------------------------------------------------------------!

    yparm =  list(0 for i in range(ne))

    RK1 = f_sys( x, y )        # ----------------------!  Compute k1

    hh = h/3.0                 # ----------------------!  Compute k2
    xparm = x + hh
    for k in range(ne):
        yparm[k] = y[k] + hh*RK1[k]

    RK2 = f_sys( xparm, yparm )

    hh = h/3.0                 # ----------------------!  Compute k3
    xparm = x + 2.0*hh
    for k in range(ne):
        yparm[k] = y[k] + hh*( -RK1[k] + 3.0*RK2[k] )

    RK3 = f_sys( xparm, yparm )

    hh = h                      # ----------------------!  Compute k4
    xparm = x + hh;
    for k in range(ne):
        yparm[k] = y[k] + hh*( RK1[k] - RK2[k] + RK3[k] )

    RK4 = f_sys( xparm, yparm )

    hh = h/8.0                  # ----------------------!  sum it
    x += h
    for k in range(ne):
        y[k] = y[k] + hh*( RK1[k] + 3.0*RK2[k] + 3.0*RK3[k] + RK4[k] )

    return x, y


#! =============================================================================!
#|
#|   Function:  RK5th
#|
#|   Purpose:   5th order Runge-Kutta-Butcher's method to solve ODE system
#|
#|   Usage:     RK5th( ne, f_sys, x, y, h )
#|
#|   Arguments:
#|     ne     - number of ODE system equation [i]
#|     f_sys  - subroutine for ODE system     [i]
#|     x      - ODE independent variable      [i/o]
#|              it is replaced to next step
#|     y      - ODE dependent variable vector [i/o]
#|     h      - step size of increment        [i]
#|  // pole   - status for the solution approach to a stiff pole [o]
#|
#|   Required:  Call system equation subroutine "f_sys" in order to input
#|              ODE system
#|
#|   Reference: 1. Curtis F. Gerald, Patrick O. Wheatley,
#|                "Applied Numerical Analysis", 4e, Addison Wesley, (1989)
#|              2. S.C.Chapra, R.P.Canale,
#|                "Numerical Methods for Engineers", 4e, (2003)
#|
#!============================================================================!
###
def RK5th( ne, f_sys, x, y, h ):

    # ----------------------------------------------------------------!
    #     Start procedure
    # ----------------------------------------------------------------!

    yparm =  list(0 for i in range(ne))

    RK1 = f_sys( x, y )      # ----------------------!  Compute k1

    hh = h/4.0               # ----------------------!  Compute k2
    xparm = x + hh
    for k in range(ne):
        yparm[k] = y[k] + hh*RK1[k]

    RK2 = f_sys( xparm, yparm )

    hh = h/8.0               # ----------------------!  Compute k3
    xparm = x + 2.0*hh;
    for k in range(ne):
        yparm[k] = y[k] + hh*( RK1[k] + RK2[k] )

    RK3 = f_sys( xparm, yparm )

    hh = h/2.0               # ----------------------!  Compute k4
    xparm = x + hh
    for k in range(ne):
        yparm[k] = y[k] + hh*( -RK2[k] + 2.0*RK3[k] )

    RK4 = f_sys( xparm, yparm )

    hh = h/16.0              # ----------------------!  Compute k5
    xparm = x + 8.0*hh
    for k in range(ne):
        yparm[k] = y[k] + hh*( 3.0*RK1[k] + 9.0*RK4[k] )
      # yparm[k] = y[k] + hh*( 3.0*RK1[k] +10.0*RK4[k] )

    RK5 = f_sys( xparm, yparm )

    hh = h/7.0               #  ----------------------!  Compute k5
    xparm = x + h
    for k in range(ne):
        yparm[k] = y[k] + hh*( -  3.0*RK1[k] +  2.0*RK2[k]  \
                               + 12.0*RK3[k] - 12.0*RK4[k]  \
                               +  8.0*RK5[k] )

    RK6 = f_sys( xparm, yparm )

    hh = h/90.0              # ----------------------!  sum it
    x += h
    for k in range(ne):
        y[k] = y[k] + hh*(   7.0*RK1[k] + 32.0*RK3[k]  \
                           + 7.0*RK6[k] + 32.0*RK5[k]  \
                           +12.0*RK4[k] )
            
    return x, y



#|============================================================================!
#|
#|   Function:  Shoot_Bisec (double precision version)                              
#|
#|   Coding:    BeH Research Lab.
#|   Compiler:  Python 3.0 over
#|
#|   Purpose:   Find unknown initial value by "shooting method"                 
#|              from given boudary value with bisection method.                 
#|                                                                              
#|   Usage:     Shoot_Bisec( ne, a, b, h, bv, eps, f_sys )
#|                                                                              
#|   Arguments:                                                                 
#|     ne     - number of ODE system equations                                  
#|     a      - starting (initial) point [i]                                    
#|     b      - ending  (boundary) point [i]                                    
#|     h      - increment of independent variable from a to b [i]               
#|     bv[i]  - array of boundary values and returned initial value [i/o]
#|     eps    - error criterion for finding exact unknown IC [i]
#|     f_sys  - name of subroutine which support the system equation [i]
#|                                                                              
#|     nc     - number of unknown I.C. (starting values)
#|     f(i)   - array of I.C. vector [i/o] : function of x "f(x)"                 
#|     z(i)   - array of unknown I.C. vector
#|              behaves as independent variable of function g(z)                
#|     g(i)   - array of   known B.C. vector [i]                                  
#|              behaves as function of z "g(z)"                                 
#|                                                                              
#|   Reference: Author                                                          
#|            - Reference Name                                                  
#|                                                                              
#!============================================================================!
###
import os

def Shoot_Bisec( ne, a, b, h, bv, eps, f_sys ):

    #----------------------------------------------------------------!
    #     Compute 2 initial point values of function g(z)
    #----------------------------------------------------------------!
    maxit = 200
    f =  list(0 for i in range(ne))
    
    f0 = f[0] = bv[0]                #  save f(0)
    f1 = f[1] = bv[1]                #  save f'(0)
    f2 = f[2] = 0.0                  #  f''(0) by arbitrary 1

    imax = (int)((b-a)/h) 
    x  = a
    z1 = f[2]
    for i in range(imax):            #  get the solution at boundary point
        x, f = RK5th( ne, f_sys, x, f, h )
        
    g1 = bv[2] - f[1]                #  f[1] = f'(inf) after compute

    f2 = f[2] = 1.0                  #  f''(0) by arbitrary 2
    x  = a
    z2 = f[2]
    for i in range(imax):            #  get the solution at boundary point
        x, f = RK5th( ne, f_sys, x, f, h )

    g2 = bv[2] - f[1]                #  f[1] = f'(inf) after compute
    
    if g1*g2 >= 0.0:  
        print( "*** No Bracket, Set I.C. again..." )
        os.system("pause")

    #----------------------------------------------------------------!
    #     Start bisection
    #----------------------------------------------------------------!

    if  g1 < 0.0:
        zz = z1;  dz = z2 - z1;
    else:
        zz = z2;  dz = z1 - z2;

    iter = 0
    while abs(dz) > eps/100.0:
        iter += 1
        if  iter > maxit:
            print( "*** Too many bisections in Root_Bisec..." )
            exit(0)
        
        dz = dz/2.0
        zm = zz + dz

        x    = a
        f[0] = f0;  f[1] = f1; f[2] = zm;   # new bisectered value
        for i in range(imax):
        #   print(f" x : {x:3.2f}    f : {f[0]: 13.9f}, {f[1]: 13.9f}, {f[2]: 13.9f} " )
            x, f = RK5th( ne, f_sys, x, f, h )

        g2 = bv[2] - f[1];
        if  g2 <= 0.0:  zz = zm
      
    print( f"      Max. epsilon       = {eps: 8.2e}" )
    print( f"      Bisction iteration = {iter: 4d}"  )


    #----------------------------------------------------------------!
    #     Restore initial values & Clear memory allocation
    #----------------------------------------------------------------!

    bv[0] = f0;  bv[1] = f1; bv[2] = zz;

    return bv   # retrurn inital value



#|=============================================================================!
#|
#|   Function:  Shoot_NR (double precision version)
#|                                                                              
#|   Coding:    BeH Research Lab.
#|   Compiler:  Python 3.0 over
#|
#|   Purpose:   Find unknown initial value by "shooting method"
#|              from given boudary value with Newton-Rhapson method.
#|
#|   Usage:     Shoot_NR( ne, a, b, h, bv, eps, f_sys )                       
#|
#|   Arguments:                                                                 
#|     ne     - number of ODE system equations [i]
#|     a      - starting (initial) point [i]
#|     b      - ending  (boundary) point [i]
#|     h      - increment of independent variable from a to b [i]
#|     bv[i]  - array of boundary value and returned initial values [i/o]
#|     eps    - error criterion for finding exact unknown IC [i]
#|     f_sys  - name of subroutine which support the system equation [i]
#|
#|     f[i]   - array of IC vector : function of x "f(x)"                 
#|     z      - unknown IC: returned as computed value
#|              behaves as independent variable of function g(z)
#|     g      - unknown BC function : behaves as function of z "g(z)"
#|
#|   Reference: Author
#|            - reference Name
#|
#!============================================================================!
###
def Shoot_NR( ne, a, b, h, bv, eps, f_sys ):

    #----------------------------------------------------------------!
    #     Compute 1st IC w/ guess value
    #----------------------------------------------------------------!
    
    f =  list(0 for i in range(ne))

    f0 = f[0] = bv[0]                #  save f(0)
    f1 = f[1] = bv[1]                #  save f'(0)
    f2 = f[2] = 0.0                  #  f''(0) by arbitrary

    imax = (int)((b-a)/h)
    x  = a
    for i in range(imax):            #  get the solution at boundary point
        x, f = RK5th( ne, f_sys, x, f, h )

    g  = error = bv[2] - f[1];       #  f[1] = f'(inf) after compute
    dgdz  = 0.1;                     #  any increment at 1st time
    z0 = f2; g0 = g;                 #  substitute to old value
    z = f2 + dgdz;                   #  set unknown initial value


    #----------------------------------------------------------------!
    #     Change initial value and check error
    #----------------------------------------------------------------!

    iter = 0;
    while abs(error) > eps:
        x = a
        f[0] = f0; f[1] = f1; f[2] = z;    #  set IC with new value: z
        for i in range(imax):
            x, f = RK5th( ne, f_sys, x, f, h )
            # x = x + h

        g = error = bv[2] - f[1]           # bv[2] = f'(inf)
        dgdz = g / (g - g0) * (z - z0)
        z0 = z;  g0 = g                    # substitute to old value
        z -= dgdz                          # Newton's iteration
        iter += 1
      
    print( f"      Max. epsilon     = {eps: 8.2e}" )
    print( f"      Newton iteration = {iter: 4d} " )


    #----------------------------------------------------------------!
    #     Restore initial values & Clear memory allocation
    #----------------------------------------------------------------!

    bv[0] = f0;  bv[1] = f1; bv[2] = z;

    return bv   # retrurn inital value



#|=============================================================================!
#|
#|   Function:  Write_Result (double precision version)
#|
#|   Purpose:   print out the results to given file pointer
#|
#|   Usage:     Write_Result( fp, method, eta, f, flag )
#|
#|   Arguments:
#|     fp      - file pointer [i]
#|     method  - name of solving method [i]
#|     eta     - independent variable [i]
#|     f[i]    - dependent data array [i]
#|     flag    - flag for header printing [i]
#|
#|============================================================================!
###
def Write_Result( fp, method, eta, f, flag ):

    if  flag == 1:
        fp.write(f"\n         The numerical solution of Blasius equation" )
        fp.write(f"\n             shooting by {method}")
        fp.write(f"\n\n\n     " )
        fp.write(f"eta      f(eta)        f'(eta)       f''(eta)         V    \n   " )
        fp.write(f"--------------------------------------------------------------\n" )
      

    fp.write("%9.4f%14.9f%14.9f%14.9f%14.9f\n"   \
              %(eta,f[0],f[1],f[2],(eta*f[1]-f[0])/2.0) )
    
    return



#!----------------------------------------------------------------------------!
#|   Main Function 
#!----------------------------------------------------------------------------!
###
import os

def main():
    
    #----------------------------------------------------------------!
    #     Set condition parameters & Open output file
    #----------------------------------------------------------------!

    ne = 3; nc = 1;          # number of system equation and unknown I.C.
    a =  0.00                # low  limit of eta
    b = 12.00                # high limit of eta
    h =  0.01                # increment of eta value
    f = [0.0, 0.0, 1.0]      # [f(0), f'(0), f''(0)]

    eps  = 0.1e-9
    iprn = 10
       
    #----------------------------------------------------------------!
    #     Choose method & Execute shooting to find initial values
    #----------------------------------------------------------------!

    method  = "Newton-Rhapson method"
    outfile = "Blasius_01.out"
#     method  = "Bisection method"
#     outfile = "Blasius_02.out"

    print(f"      Shooting by {method}" )
    print(f"      output file : \'{outfile}\'\n" )

    f = Shoot_NR( ne, a, b, h, f, eps, BlasiusEqn )    
#     f = Shoot_Bisec( ne, a, b, h, f, eps, BlasiusEqn )    
   
    #----------------------------------------------------------------!
    #     Save initial value & Write result
    #----------------------------------------------------------------!
    
    fp = open( outfile, 'w' )    
    eta = a
    f0 = f[0]; f1 = f[1]; f2 = f[2];  # f(0), f'(0) and f''(0) after computing

    Write_Result( fp, method, eta, f, flag=1 )  # print header


    # print numerical result by eta increament
    imax = int((b-a)/h)
    for i in range(imax):

        eta, f = RK5th( ne, BlasiusEqn, eta, f, h )

        if (i+1)%iprn == 0:
            Write_Result( fp, method, eta, f, flag=0 )

    fp.close()

    print( f"\n      delta eta = {h: 7.3f}"    )
    print( f"      max. eta  = {eta : 16.12f}\n" )
    print( f"      Calculated I.C's : [{f0}, {f1}, {f2}]")
    print( f"      f(inf)   = {f[0]: 16.12f}"  )
    print( f"      f'(inf)  = {f[1]: 16.12f}")
    print( f"      f''(inf) = {f[2]: 16.12f}\n")

    print( f">>> Execution is completed..." )
    os.system( "pause" )
    return

    
if __name__ == '__main__':
    main()





