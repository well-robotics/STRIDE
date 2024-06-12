/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:39 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t2393;
  double t2431;
  double t2318;
  double t2358;
  double t2310;
  double t2350;
  double t2395;
  double t2402;
  double t2417;
  double t2418;
  double t2419;
  double t2441;
  double t2442;
  double t2443;
  double t2447;
  double t2448;
  double t2449;
  double t2450;
  double t2451;
  double t2452;
  double t2453;
  double t2434;
  double t2435;
  double t2436;
  double t2437;
  double t2438;
  double t2439;
  double t2455;
  double t2460;
  double t2461;
  double t2462;
  double t2463;
  double t2454;
  double t2456;
  double t2464;
  double t2319;
  double t2376;
  double t2377;
  double t2378;
  double t2444;
  double t2476;
  double t2477;
  double t2478;
  double t2500;
  double t2501;
  double t2502;
  double t2492;
  double t2493;
  double t2494;
  double t2496;
  double t2497;
  double t2498;
  double t2517;
  double t2518;
  double t2519;
  double t2521;
  double t2522;
  double t2523;
  double t2440;
  double t2445;
  double t2475;
  double t2479;
  double t2480;
  double t2481;
  double t2484;
  double t2485;
  double t2486;
  double t2487;
  double t2488;
  double t2489;
  double t2499;
  double t2503;
  double t2547;
  double t2548;
  double t2508;
  double t2550;
  double t2551;
  double t2510;
  t2393 = Cos(var1[6]);
  t2431 = Sin(var1[6]);
  t2318 = Sin(var1[2]);
  t2358 = Sin(var1[5]);
  t2310 = Cos(var1[5]);
  t2350 = Cos(var1[2]);
  t2395 = -1.*t2393;
  t2402 = 1. + t2395;
  t2417 = 0.5*t2402;
  t2418 = 0.671885*t2393;
  t2419 = t2417 + t2418;
  t2441 = t2310*t2393;
  t2442 = -1.*t2358*t2431;
  t2443 = t2441 + t2442;
  t2447 = t2419*t2393;
  t2448 = Power(t2431,2);
  t2449 = 0.171885*t2448;
  t2450 = t2447 + t2449;
  t2451 = -1.*t2393*t2358;
  t2452 = -1.*t2310*t2431;
  t2453 = t2451 + t2452;
  t2434 = t2419*t2431;
  t2435 = -0.171885*t2393*t2431;
  t2436 = t2434 + t2435;
  t2437 = t2393*t2358;
  t2438 = t2310*t2431;
  t2439 = t2437 + t2438;
  t2455 = -1.*t2318*t2443;
  t2460 = -1.*t2350*t2310;
  t2461 = t2318*t2358;
  t2462 = t2460 + t2461;
  t2463 = -0.51185934*t2462;
  t2454 = t2350*t2453;
  t2456 = t2454 + t2455;
  t2464 = -1.*t2318*t2453;
  t2319 = -1.*t2310*t2318;
  t2376 = -1.*t2350*t2358;
  t2377 = t2319 + t2376;
  t2378 = -0.51185934*t2377;
  t2444 = t2350*t2443;
  t2476 = -1.*t2310*t2393;
  t2477 = t2358*t2431;
  t2478 = t2476 + t2477;
  t2500 = t2310*t2419;
  t2501 = -0.171885*t2358*t2431;
  t2502 = t2500 + t2501;
  t2492 = -1.*t2419*t2358;
  t2493 = -0.171885*t2310*t2431;
  t2494 = t2492 + t2493;
  t2496 = t2419*t2358;
  t2497 = 0.171885*t2310*t2431;
  t2498 = t2496 + t2497;
  t2517 = -1.*t2419*t2431;
  t2518 = 0.171885*t2393*t2431;
  t2519 = t2517 + t2518;
  t2521 = Power(t2393,2);
  t2522 = -0.171885*t2521;
  t2523 = t2447 + t2522;
  t2440 = -1.*t2318*t2439;
  t2445 = t2440 + t2444;
  t2475 = -0.85216*t2436*t2456;
  t2479 = t2350*t2478;
  t2480 = t2464 + t2479;
  t2481 = -0.85216*t2450*t2480;
  t2484 = t2318*t2453;
  t2485 = t2484 + t2444;
  t2486 = -0.85216*t2436*t2485;
  t2487 = t2318*t2478;
  t2488 = t2454 + t2487;
  t2489 = -0.85216*t2450*t2488;
  t2499 = t2498*t2443;
  t2503 = t2453*t2502;
  t2547 = -0.171885*t2393*t2358;
  t2548 = t2547 + t2493;
  t2508 = -1.*t2453*t2498;
  t2550 = 0.171885*t2310*t2393;
  t2551 = t2550 + t2501;
  t2510 = -1.*t2502*t2478;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(t2378 - 0.85216*t2436*t2445 - 0.85216*t2450*t2456)*var2[0] - 0.5*(-0.85216*t2436*(-1.*t2350*t2439 + t2455) + t2463 - 0.85216*t2450*(-1.*t2350*t2443 + t2464))*var2[1])*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*(t2378 + t2486 + t2489)*var2[0] - 0.5*(t2463 + t2475 + t2481)*var2[1] - 0.5*(-0.85216*t2450*(t2443*t2494 + t2499 + t2439*t2502 + t2503) - 0.85216*t2436*(-1.*t2453*t2494 - 1.*t2443*t2502 + t2508 + t2510))*var2[2])*var2[5];
  p_output1[6]=var2[5]*(-0.5*(t2486 + t2489 - 0.85216*t2485*t2519 - 0.85216*(t2350*t2439 + t2318*t2443)*t2523)*var2[0] - 0.5*(t2475 + t2481 - 0.85216*t2456*t2519 - 0.85216*t2445*t2523)*var2[1] - 0.5*(-0.85216*(t2439*t2498 + t2443*t2502)*t2519 - 0.85216*(-1.*t2443*t2498 - 1.*t2453*t2502)*t2523 - 0.85216*t2450*(t2499 + t2503 + t2443*t2548 + t2439*t2551) - 0.85216*t2436*(t2508 + t2510 - 1.*t2453*t2548 - 1.*t2443*t2551))*var2[2] - 0.5*(-1.70432*t2450*t2519 - 1.70432*t2436*t2523)*var2[5] + 0.0732367608*t2519*var2[6]);
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

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
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
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce3_vec6_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
