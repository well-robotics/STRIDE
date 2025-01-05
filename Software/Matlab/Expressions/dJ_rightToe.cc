/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:48 GMT-05:00
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
  double t2190;
  double t2276;
  double t2295;
  double t2313;
  double t282;
  double t2336;
  double t2357;
  double t2361;
  double t2364;
  double t2377;
  double t2350;
  double t2351;
  double t2355;
  double t2398;
  double t2399;
  double t2404;
  double t2325;
  double t2365;
  double t2384;
  double t2390;
  double t2406;
  double t2409;
  double t2413;
  double t2464;
  double t2469;
  double t2470;
  double t2475;
  double t2476;
  double t2479;
  double t2483;
  double t2489;
  double t2495;
  double t2497;
  double t2499;
  double t2501;
  double t2502;
  double t2505;
  double t2435;
  double t2437;
  double t2438;
  double t2300;
  double t2329;
  double t2341;
  double t2342;
  double t2344;
  double t2424;
  double t2488;
  double t2445;
  double t2449;
  double t2451;
  double t2457;
  double t2458;
  double t2506;
  double t2538;
  double t2539;
  double t2395;
  double t2420;
  double t2428;
  double t2431;
  double t2432;
  double t2439;
  double t2450;
  double t2463;
  double t2474;
  double t2480;
  double t2484;
  double t2485;
  double t2486;
  double t2490;
  double t2492;
  double t2493;
  double t2572;
  double t2573;
  double t2574;
  double t2536;
  double t2537;
  double t2540;
  double t2543;
  double t2544;
  double t2545;
  double t2549;
  double t2550;
  double t2551;
  double t2552;
  double t2553;
  double t2555;
  double t2556;
  double t2557;
  double t2560;
  double t2563;
  double t2564;
  double t2565;
  double t2566;
  double t2567;
  double t2568;
  double t2569;
  double t2498;
  double t2503;
  double t2507;
  double t2508;
  double t2511;
  double t2512;
  double t2514;
  double t2592;
  double t2593;
  double t2520;
  double t2522;
  double t2523;
  double t2526;
  double t2527;
  double t2530;
  double t2531;
  double t2589;
  double t2591;
  double t2594;
  double t2595;
  double t2596;
  double t2602;
  double t2603;
  double t2604;
  t2190 = Cos(var1[3]);
  t2276 = -1.*t2190;
  t2295 = 1. + t2276;
  t2313 = Sin(var1[3]);
  t282 = Cos(var1[2]);
  t2336 = Sin(var1[2]);
  t2357 = Cos(var1[4]);
  t2361 = -1.*t2357;
  t2364 = 1. + t2361;
  t2377 = Sin(var1[4]);
  t2350 = -1.*t282*t2190;
  t2351 = -1.*t2336*t2313;
  t2355 = t2350 + t2351;
  t2398 = -1.*t2190*t2336;
  t2399 = t282*t2313;
  t2404 = t2398 + t2399;
  t2325 = -0.0695*t2313;
  t2365 = -0.0265*t2364;
  t2384 = -0.2375*t2377;
  t2390 = t2365 + t2384;
  t2406 = -0.2375*t2364;
  t2409 = 0.0265*t2377;
  t2413 = t2406 + t2409;
  t2464 = t282*t2190;
  t2469 = t2336*t2313;
  t2470 = t2464 + t2469;
  t2475 = t2190*t2336;
  t2476 = -1.*t282*t2313;
  t2479 = t2475 + t2476;
  t2483 = t2357*t2470;
  t2489 = -1.*t2470*t2377;
  t2495 = 0.0265*t2357;
  t2497 = t2495 + t2384;
  t2499 = -0.2375*t2357;
  t2501 = -0.0265*t2377;
  t2502 = t2499 + t2501;
  t2505 = -1.*t2404*t2377;
  t2435 = t2357*t2404;
  t2437 = -1.*t2355*t2377;
  t2438 = t2435 + t2437;
  t2300 = -0.0265*t2295;
  t2329 = t2300 + t2325;
  t2341 = -0.0695*t2295;
  t2342 = 0.0265*t2313;
  t2344 = t2341 + t2342;
  t2424 = t2357*t2355;
  t2488 = t2357*t2479;
  t2445 = 0.0265*t2190;
  t2449 = t2445 + t2325;
  t2451 = -0.0695*t2190;
  t2457 = -0.0265*t2313;
  t2458 = t2451 + t2457;
  t2506 = t2483 + t2505;
  t2538 = -1.*t2479*t2377;
  t2539 = t2424 + t2538;
  t2395 = t2355*t2390;
  t2420 = t2404*t2413;
  t2428 = t2404*t2377;
  t2431 = t2424 + t2428;
  t2432 = -0.0265*t2431;
  t2439 = -0.0115*t2438;
  t2450 = t282*t2449;
  t2463 = -1.*t2336*t2458;
  t2474 = t2470*t2390;
  t2480 = t2479*t2413;
  t2484 = t2479*t2377;
  t2485 = t2483 + t2484;
  t2486 = -0.0265*t2485;
  t2490 = t2488 + t2489;
  t2492 = -0.0115*t2490;
  t2493 = t2450 + t2463 + t2474 + t2480 + t2486 + t2492;
  t2572 = -0.0265*t2190;
  t2573 = 0.0695*t2313;
  t2574 = t2572 + t2573;
  t2536 = t2479*t2390;
  t2537 = t2355*t2413;
  t2540 = -0.0115*t2539;
  t2543 = t2355*t2377;
  t2544 = t2488 + t2543;
  t2545 = -0.0265*t2544;
  t2549 = -1.*t2336*t2449;
  t2550 = -1.*t282*t2458;
  t2551 = t2404*t2390;
  t2552 = t2470*t2413;
  t2553 = -0.0115*t2506;
  t2555 = t2470*t2377;
  t2556 = t2435 + t2555;
  t2557 = -0.0265*t2556;
  t2560 = t2549 + t2550 + t2551 + t2552 + t2553 + t2557;
  t2563 = t2355*t2497;
  t2564 = t2479*t2502;
  t2565 = -0.0265*t2539;
  t2566 = -1.*t2357*t2479;
  t2567 = t2566 + t2437;
  t2568 = -0.0115*t2567;
  t2569 = t2563 + t2564 + t2565 + t2568;
  t2498 = t2470*t2497;
  t2503 = t2404*t2502;
  t2507 = -0.0265*t2506;
  t2508 = -1.*t2357*t2404;
  t2511 = t2508 + t2489;
  t2512 = -0.0115*t2511;
  t2514 = t2498 + t2503 + t2507 + t2512;
  t2592 = -1.*t2357*t2470;
  t2593 = t2592 + t2538;
  t2520 = t2404*t2497;
  t2522 = t2355*t2502;
  t2523 = -1.*t2357*t2355;
  t2526 = t2523 + t2505;
  t2527 = -0.0115*t2526;
  t2530 = -0.0265*t2438;
  t2531 = t2520 + t2522 + t2527 + t2530;
  t2589 = t2479*t2497;
  t2591 = t2470*t2502;
  t2594 = -0.0115*t2593;
  t2595 = -0.0265*t2490;
  t2596 = t2589 + t2591 + t2594 + t2595;
  t2602 = -0.0265*t2357;
  t2603 = 0.2375*t2377;
  t2604 = t2602 + t2603;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=(-1.*t2336*t2344 + t2395 + t2420 + t2432 + t2439 - 1.*t2329*t282)*var1[9] + t2493*var1[10] + t2514*var1[11];
  p_output1[5]=(t2329*t2336 + t2536 + t2537 + t2540 + t2545 - 1.*t2344*t282)*var1[9] + t2560*var1[10] + t2531*var1[11];
  p_output1[6]=t2493*var1[9] + (t2395 + t2420 + t2432 + t2439 + t2336*t2458 + t2574*t282)*var1[10] + t2569*var1[11];
  p_output1[7]=t2560*var1[9] + (t2536 + t2537 + t2540 + t2545 - 1.*t2336*t2574 + t2458*t282)*var1[10] + t2596*var1[11];
  p_output1[8]=t2514*var1[9] + t2569*var1[10] + (t2564 - 0.0115*(t2555 + t2566) - 0.0265*t2593 + t2470*t2604)*var1[11];
  p_output1[9]=t2531*var1[9] + t2596*var1[10] + (-0.0265*t2511 + t2591 - 0.0115*(t2428 + t2592) + t2404*t2604)*var1[11];
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
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
    ( !(mrows == 14 && ncols == 1) && 
      !(mrows == 1 && ncols == 14))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "dJ_rightToe.hh"

namespace SymFunction
{

void dJ_rightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
