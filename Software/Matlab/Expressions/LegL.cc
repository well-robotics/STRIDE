/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:01 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
#include<math.h>
/**
 * Copied from Wolfram Mathematica C Definitions file mdefs.hpp
 * Changed marcos to inline functions (Eric Cousineau)
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

/* update ArcTan function to use atan2 instead. */
inline double ArcTan(double x, double y) { return atan2(y,x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;

inline double Sec(double x) { return 1/cos(x); }
inline double Csc(double x) { return 1/sin(x); }

#endif

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t2325;
  double t2334;
  double t2341;
  double t2346;
  double t2300;
  double t2399;
  double t2438;
  double t2440;
  double t2444;
  double t2450;
  double t2428;
  double t2431;
  double t2435;
  double t2469;
  double t2474;
  double t2475;
  double t2342;
  double t2365;
  double t2384;
  double t2406;
  double t2409;
  double t2413;
  double t2445;
  double t2451;
  double t2457;
  double t2480;
  double t2483;
  double t2484;
  double t2486;
  double t2506;
  double t2507;
  double t2511;
  double t2544;
  double t2546;
  double t2548;
  double t2552;
  double t2568;
  double t2569;
  double t2570;
  double t2572;
  double t2564;
  double t2565;
  double t2567;
  double t2581;
  double t2582;
  double t2583;
  double t2551;
  double t2553;
  double t2555;
  double t2557;
  double t2561;
  double t2562;
  double t2571;
  double t2573;
  double t2579;
  double t2584;
  double t2585;
  double t2587;
  double t2589;
  double t2604;
  double t2605;
  double t2606;
  t2325 = Cos(var1[3]);
  t2334 = -1.*t2325;
  t2341 = 1. + t2334;
  t2346 = Sin(var1[3]);
  t2300 = Cos(var1[2]);
  t2399 = Sin(var1[2]);
  t2438 = Cos(var1[4]);
  t2440 = -1.*t2438;
  t2444 = 1. + t2440;
  t2450 = Sin(var1[4]);
  t2428 = t2300*t2325;
  t2431 = t2399*t2346;
  t2435 = t2428 + t2431;
  t2469 = t2325*t2399;
  t2474 = -1.*t2300*t2346;
  t2475 = t2469 + t2474;
  t2342 = -0.0265*t2341;
  t2365 = -0.0695*t2346;
  t2384 = t2342 + t2365;
  t2406 = -0.0695*t2341;
  t2409 = 0.0265*t2346;
  t2413 = t2406 + t2409;
  t2445 = -0.0265*t2444;
  t2451 = -0.2375*t2450;
  t2457 = t2445 + t2451;
  t2480 = -0.2375*t2444;
  t2483 = 0.0265*t2450;
  t2484 = t2480 + t2483;
  t2486 = t2438*t2435;
  t2506 = -1.*t2325*t2399;
  t2507 = t2300*t2346;
  t2511 = t2506 + t2507;
  t2544 = Cos(var1[5]);
  t2546 = -1.*t2544;
  t2548 = 1. + t2546;
  t2552 = Sin(var1[5]);
  t2568 = Cos(var1[6]);
  t2569 = -1.*t2568;
  t2570 = 1. + t2569;
  t2572 = Sin(var1[6]);
  t2564 = t2300*t2544;
  t2565 = -1.*t2399*t2552;
  t2567 = t2564 + t2565;
  t2581 = -1.*t2544*t2399;
  t2582 = -1.*t2300*t2552;
  t2583 = t2581 + t2582;
  t2551 = -0.0695*t2548;
  t2553 = -0.0265*t2552;
  t2555 = t2551 + t2553;
  t2557 = -0.0265*t2548;
  t2561 = 0.0695*t2552;
  t2562 = t2557 + t2561;
  t2571 = -0.2375*t2570;
  t2573 = -0.0265*t2572;
  t2579 = t2571 + t2573;
  t2584 = -0.0265*t2570;
  t2585 = 0.2375*t2572;
  t2587 = t2584 + t2585;
  t2589 = t2568*t2567;
  t2604 = t2544*t2399;
  t2605 = t2300*t2552;
  t2606 = t2604 + t2605;
  p_output1[0]=Sqrt(0.00085849 + Power(-1.*t2300*t2384 - 1.*t2399*t2413 - 1.*t2435*t2457 + 0.0115*(-1.*t2435*t2450 + t2438*t2475) - 1.*t2475*t2484 + 0.0265*(t2450*t2475 + t2486),2) + Power(t2384*t2399 - 1.*t2300*t2413 - 1.*t2435*t2484 - 1.*t2457*t2511 + 0.0265*(t2435*t2450 + t2438*t2511) + 0.0115*(t2486 - 1.*t2450*t2511),2));
  p_output1[1]=Sqrt(0.00085849 + Power(-1.*t2300*t2555 + t2399*t2562 - 1.*t2567*t2579 + 0.0265*(-1.*t2567*t2572 + t2568*t2583) - 1.*t2583*t2587 + 0.0115*(t2572*t2583 + t2589),2) + Power(-1.*t2399*t2555 - 1.*t2300*t2562 - 1.*t2567*t2587 - 1.*t2579*t2606 + 0.0115*(t2567*t2572 + t2568*t2606) + 0.0265*(t2589 - 1.*t2572*t2606),2));
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "LegL.hh"

namespace SymFunction
{

void LegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
