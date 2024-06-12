/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:51 GMT-05:00
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

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1)
{
  double t2272;
  double t2279;
  double t2285;
  double t2293;
  double t2256;
  double t2347;
  double t2382;
  double t2387;
  double t2388;
  double t2394;
  double t2374;
  double t2375;
  double t2380;
  double t2413;
  double t2418;
  double t2419;
  double t2287;
  double t2320;
  double t2333;
  double t2352;
  double t2356;
  double t2363;
  double t2392;
  double t2400;
  double t2401;
  double t2426;
  double t2427;
  double t2428;
  double t2431;
  double t2450;
  double t2451;
  double t2455;
  double t2488;
  double t2491;
  double t2492;
  double t2496;
  double t2512;
  double t2513;
  double t2514;
  double t2516;
  double t2508;
  double t2509;
  double t2511;
  double t2525;
  double t2526;
  double t2527;
  double t2495;
  double t2498;
  double t2499;
  double t2503;
  double t2505;
  double t2506;
  double t2515;
  double t2517;
  double t2523;
  double t2528;
  double t2530;
  double t2531;
  double t2534;
  double t2548;
  double t2549;
  double t2550;
  t2272 = Cos(var1[3]);
  t2279 = -1.*t2272;
  t2285 = 1. + t2279;
  t2293 = Sin(var1[3]);
  t2256 = Cos(var1[2]);
  t2347 = Sin(var1[2]);
  t2382 = Cos(var1[4]);
  t2387 = -1.*t2382;
  t2388 = 1. + t2387;
  t2394 = Sin(var1[4]);
  t2374 = t2256*t2272;
  t2375 = t2347*t2293;
  t2380 = t2374 + t2375;
  t2413 = t2272*t2347;
  t2418 = -1.*t2256*t2293;
  t2419 = t2413 + t2418;
  t2287 = -0.0265*t2285;
  t2320 = -0.0695*t2293;
  t2333 = t2287 + t2320;
  t2352 = -0.0695*t2285;
  t2356 = 0.0265*t2293;
  t2363 = t2352 + t2356;
  t2392 = -0.0265*t2388;
  t2400 = -0.2375*t2394;
  t2401 = t2392 + t2400;
  t2426 = -0.2375*t2388;
  t2427 = 0.0265*t2394;
  t2428 = t2426 + t2427;
  t2431 = t2382*t2380;
  t2450 = -1.*t2272*t2347;
  t2451 = t2256*t2293;
  t2455 = t2450 + t2451;
  t2488 = Cos(var1[5]);
  t2491 = -1.*t2488;
  t2492 = 1. + t2491;
  t2496 = Sin(var1[5]);
  t2512 = Cos(var1[6]);
  t2513 = -1.*t2512;
  t2514 = 1. + t2513;
  t2516 = Sin(var1[6]);
  t2508 = t2256*t2488;
  t2509 = -1.*t2347*t2496;
  t2511 = t2508 + t2509;
  t2525 = -1.*t2488*t2347;
  t2526 = -1.*t2256*t2496;
  t2527 = t2525 + t2526;
  t2495 = -0.0695*t2492;
  t2498 = -0.0265*t2496;
  t2499 = t2495 + t2498;
  t2503 = -0.0265*t2492;
  t2505 = 0.0695*t2496;
  t2506 = t2503 + t2505;
  t2515 = -0.2375*t2514;
  t2517 = -0.0265*t2516;
  t2523 = t2515 + t2517;
  t2528 = -0.0265*t2514;
  t2530 = 0.2375*t2516;
  t2531 = t2528 + t2530;
  t2534 = t2512*t2511;
  t2548 = t2488*t2347;
  t2549 = t2256*t2496;
  t2550 = t2548 + t2549;
  p_output1[0]=Sqrt(0.00085849 + Power(-1.*t2256*t2333 - 1.*t2347*t2363 - 1.*t2380*t2401 - 0.0325*(-1.*t2380*t2394 + t2382*t2419) - 1.*t2419*t2428 + 0.0265*(t2394*t2419 + t2431),2) + Power(t2333*t2347 - 1.*t2256*t2363 - 1.*t2380*t2428 - 1.*t2401*t2455 + 0.0265*(t2380*t2394 + t2382*t2455) - 0.0325*(t2431 - 1.*t2394*t2455),2));
  p_output1[1]=Sqrt(0.00085849 + Power(-1.*t2256*t2499 + t2347*t2506 - 1.*t2511*t2523 + 0.0265*(-1.*t2511*t2516 + t2512*t2527) - 1.*t2527*t2531 - 0.0325*(t2516*t2527 + t2534),2) + Power(-1.*t2347*t2499 - 1.*t2256*t2506 - 1.*t2511*t2531 - 1.*t2523*t2550 - 0.0325*(t2511*t2516 + t2512*t2550) + 0.0265*(t2534 - 1.*t2516*t2550),2));
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
