/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:46 GMT-05:00
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
  double t2219;
  double t2237;
  double t2243;
  double t2268;
  double t2133;
  double t2284;
  double t2304;
  double t2307;
  double t2308;
  double t2327;
  double t2294;
  double t2298;
  double t2300;
  double t2342;
  double t2347;
  double t2349;
  double t2272;
  double t2320;
  double t2333;
  double t2338;
  double t2352;
  double t2356;
  double t2363;
  double t2412;
  double t2413;
  double t2417;
  double t2419;
  double t2422;
  double t2423;
  double t2427;
  double t2433;
  double t2440;
  double t2441;
  double t2444;
  double t2445;
  double t2446;
  double t2449;
  double t2380;
  double t2381;
  double t2382;
  double t2256;
  double t2273;
  double t2285;
  double t2287;
  double t2289;
  double t2371;
  double t2432;
  double t2392;
  double t2393;
  double t2400;
  double t2401;
  double t2406;
  double t2450;
  double t2482;
  double t2483;
  double t2341;
  double t2367;
  double t2374;
  double t2375;
  double t2378;
  double t2383;
  double t2394;
  double t2407;
  double t2418;
  double t2426;
  double t2428;
  double t2429;
  double t2431;
  double t2435;
  double t2436;
  double t2437;
  double t2516;
  double t2517;
  double t2522;
  double t2480;
  double t2481;
  double t2486;
  double t2487;
  double t2488;
  double t2489;
  double t2493;
  double t2494;
  double t2495;
  double t2496;
  double t2498;
  double t2499;
  double t2500;
  double t2503;
  double t2504;
  double t2507;
  double t2508;
  double t2509;
  double t2510;
  double t2511;
  double t2512;
  double t2513;
  double t2442;
  double t2448;
  double t2451;
  double t2454;
  double t2455;
  double t2457;
  double t2458;
  double t2536;
  double t2537;
  double t2465;
  double t2466;
  double t2469;
  double t2470;
  double t2473;
  double t2474;
  double t2475;
  double t2534;
  double t2535;
  double t2538;
  double t2539;
  double t2540;
  double t2546;
  double t2547;
  double t2548;
  t2219 = Cos(var1[3]);
  t2237 = -1.*t2219;
  t2243 = 1. + t2237;
  t2268 = Sin(var1[3]);
  t2133 = Cos(var1[2]);
  t2284 = Sin(var1[2]);
  t2304 = Cos(var1[4]);
  t2307 = -1.*t2304;
  t2308 = 1. + t2307;
  t2327 = Sin(var1[4]);
  t2294 = -1.*t2133*t2219;
  t2298 = -1.*t2284*t2268;
  t2300 = t2294 + t2298;
  t2342 = -1.*t2219*t2284;
  t2347 = t2133*t2268;
  t2349 = t2342 + t2347;
  t2272 = -0.0695*t2268;
  t2320 = -0.0265*t2308;
  t2333 = -0.2375*t2327;
  t2338 = t2320 + t2333;
  t2352 = -0.2375*t2308;
  t2356 = 0.0265*t2327;
  t2363 = t2352 + t2356;
  t2412 = t2133*t2219;
  t2413 = t2284*t2268;
  t2417 = t2412 + t2413;
  t2419 = t2219*t2284;
  t2422 = -1.*t2133*t2268;
  t2423 = t2419 + t2422;
  t2427 = t2304*t2417;
  t2433 = -1.*t2417*t2327;
  t2440 = 0.0265*t2304;
  t2441 = t2440 + t2333;
  t2444 = -0.2375*t2304;
  t2445 = -0.0265*t2327;
  t2446 = t2444 + t2445;
  t2449 = -1.*t2349*t2327;
  t2380 = t2304*t2349;
  t2381 = -1.*t2300*t2327;
  t2382 = t2380 + t2381;
  t2256 = -0.0265*t2243;
  t2273 = t2256 + t2272;
  t2285 = -0.0695*t2243;
  t2287 = 0.0265*t2268;
  t2289 = t2285 + t2287;
  t2371 = t2304*t2300;
  t2432 = t2304*t2423;
  t2392 = 0.0265*t2219;
  t2393 = t2392 + t2272;
  t2400 = -0.0695*t2219;
  t2401 = -0.0265*t2268;
  t2406 = t2400 + t2401;
  t2450 = t2427 + t2449;
  t2482 = -1.*t2423*t2327;
  t2483 = t2371 + t2482;
  t2341 = t2300*t2338;
  t2367 = t2349*t2363;
  t2374 = t2349*t2327;
  t2375 = t2371 + t2374;
  t2378 = -0.0265*t2375;
  t2383 = 0.0325*t2382;
  t2394 = t2133*t2393;
  t2407 = -1.*t2284*t2406;
  t2418 = t2417*t2338;
  t2426 = t2423*t2363;
  t2428 = t2423*t2327;
  t2429 = t2427 + t2428;
  t2431 = -0.0265*t2429;
  t2435 = t2432 + t2433;
  t2436 = 0.0325*t2435;
  t2437 = t2394 + t2407 + t2418 + t2426 + t2431 + t2436;
  t2516 = -0.0265*t2219;
  t2517 = 0.0695*t2268;
  t2522 = t2516 + t2517;
  t2480 = t2423*t2338;
  t2481 = t2300*t2363;
  t2486 = 0.0325*t2483;
  t2487 = t2300*t2327;
  t2488 = t2432 + t2487;
  t2489 = -0.0265*t2488;
  t2493 = -1.*t2284*t2393;
  t2494 = -1.*t2133*t2406;
  t2495 = t2349*t2338;
  t2496 = t2417*t2363;
  t2498 = 0.0325*t2450;
  t2499 = t2417*t2327;
  t2500 = t2380 + t2499;
  t2503 = -0.0265*t2500;
  t2504 = t2493 + t2494 + t2495 + t2496 + t2498 + t2503;
  t2507 = t2300*t2441;
  t2508 = t2423*t2446;
  t2509 = -0.0265*t2483;
  t2510 = -1.*t2304*t2423;
  t2511 = t2510 + t2381;
  t2512 = 0.0325*t2511;
  t2513 = t2507 + t2508 + t2509 + t2512;
  t2442 = t2417*t2441;
  t2448 = t2349*t2446;
  t2451 = -0.0265*t2450;
  t2454 = -1.*t2304*t2349;
  t2455 = t2454 + t2433;
  t2457 = 0.0325*t2455;
  t2458 = t2442 + t2448 + t2451 + t2457;
  t2536 = -1.*t2304*t2417;
  t2537 = t2536 + t2482;
  t2465 = t2349*t2441;
  t2466 = t2300*t2446;
  t2469 = -1.*t2304*t2300;
  t2470 = t2469 + t2449;
  t2473 = 0.0325*t2470;
  t2474 = -0.0265*t2382;
  t2475 = t2465 + t2466 + t2473 + t2474;
  t2534 = t2423*t2441;
  t2535 = t2417*t2446;
  t2538 = 0.0325*t2537;
  t2539 = -0.0265*t2435;
  t2540 = t2534 + t2535 + t2538 + t2539;
  t2546 = -0.0265*t2304;
  t2547 = 0.2375*t2327;
  t2548 = t2546 + t2547;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=(-1.*t2133*t2273 - 1.*t2284*t2289 + t2341 + t2367 + t2378 + t2383)*var1[9] + t2437*var1[10] + t2458*var1[11];
  p_output1[5]=(t2273*t2284 - 1.*t2133*t2289 + t2480 + t2481 + t2486 + t2489)*var1[9] + t2504*var1[10] + t2475*var1[11];
  p_output1[6]=t2437*var1[9] + (t2341 + t2367 + t2378 + t2383 + t2284*t2406 + t2133*t2522)*var1[10] + t2513*var1[11];
  p_output1[7]=t2504*var1[9] + (t2133*t2406 + t2480 + t2481 + t2486 + t2489 - 1.*t2284*t2522)*var1[10] + t2540*var1[11];
  p_output1[8]=t2458*var1[9] + t2513*var1[10] + (t2508 + 0.0325*(t2499 + t2510) - 0.0265*t2537 + t2417*t2548)*var1[11];
  p_output1[9]=t2475*var1[9] + t2540*var1[10] + (-0.0265*t2455 + t2535 + 0.0325*(t2374 + t2536) + t2349*t2548)*var1[11];
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
