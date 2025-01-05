/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:45 GMT-05:00
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
  double t2184;
  double t2191;
  double t2141;
  double t2222;
  double t2261;
  double t2293;
  double t2190;
  double t2223;
  double t2260;
  double t2314;
  double t2322;
  double t2323;
  double t2329;
  double t2373;
  double t2376;
  double t2389;
  double t2390;
  double t2393;
  double t2338;
  double t2344;
  double t2365;
  double t2421;
  double t2423;
  double t2377;
  double t2384;
  double t2387;
  double t2395;
  double t2398;
  double t2399;
  double t2402;
  double t2404;
  double t2406;
  double t2408;
  double t2409;
  double t2413;
  double t2415;
  double t2276;
  double t2295;
  double t2300;
  double t2325;
  double t2334;
  double t2355;
  double t2357;
  double t2360;
  double t2364;
  double t2366;
  double t2457;
  double t2458;
  double t2461;
  double t2447;
  double t2424;
  double t2427;
  double t2431;
  double t2432;
  double t2434;
  double t2456;
  double t2462;
  double t2463;
  double t2464;
  double t2468;
  double t2469;
  double t2470;
  double t2472;
  double t2313;
  double t2336;
  double t2337;
  double t2341;
  double t2342;
  double t2345;
  double t2346;
  double t2350;
  double t2351;
  double t2352;
  double t2361;
  double t2372;
  double t2416;
  double t2440;
  double t2444;
  double t2445;
  double t2446;
  double t2449;
  double t2450;
  double t2451;
  double t2452;
  double t2453;
  double t2454;
  double t2473;
  double t2485;
  double t2486;
  double t2488;
  double t2503;
  double t2505;
  double t2506;
  t2184 = Cos(var1[5]);
  t2191 = Sin(var1[2]);
  t2141 = Cos(var1[2]);
  t2222 = Sin(var1[5]);
  t2261 = Cos(var1[6]);
  t2293 = Sin(var1[6]);
  t2190 = t2141*t2184;
  t2223 = -1.*t2191*t2222;
  t2260 = t2190 + t2223;
  t2314 = -1.*t2184*t2191;
  t2322 = -1.*t2141*t2222;
  t2323 = t2314 + t2322;
  t2329 = -0.0265*t2293;
  t2373 = -1.*t2261;
  t2376 = 1. + t2373;
  t2389 = -1.*t2141*t2184;
  t2390 = t2191*t2222;
  t2393 = t2389 + t2390;
  t2338 = -1.*t2323*t2293;
  t2344 = t2261*t2323;
  t2365 = -0.0265*t2222;
  t2421 = -1.*t2184;
  t2423 = 1. + t2421;
  t2377 = -0.2375*t2376;
  t2384 = t2377 + t2329;
  t2387 = t2323*t2384;
  t2395 = -0.0265*t2376;
  t2398 = 0.2375*t2293;
  t2399 = t2395 + t2398;
  t2402 = t2393*t2399;
  t2404 = t2261*t2393;
  t2406 = t2404 + t2338;
  t2408 = -0.0265*t2406;
  t2409 = t2393*t2293;
  t2413 = t2344 + t2409;
  t2415 = -0.0115*t2413;
  t2276 = -0.0265*t2261;
  t2295 = -0.2375*t2293;
  t2300 = t2276 + t2295;
  t2325 = 0.2375*t2261;
  t2334 = t2325 + t2329;
  t2355 = -0.0265*t2184;
  t2357 = -0.0695*t2222;
  t2360 = t2355 + t2357;
  t2364 = 0.0695*t2184;
  t2366 = t2364 + t2365;
  t2457 = t2184*t2191;
  t2458 = t2141*t2222;
  t2461 = t2457 + t2458;
  t2447 = -1.*t2393*t2293;
  t2424 = -0.0695*t2423;
  t2427 = t2424 + t2365;
  t2431 = -0.0265*t2423;
  t2432 = 0.0695*t2222;
  t2434 = t2431 + t2432;
  t2456 = t2393*t2384;
  t2462 = t2461*t2399;
  t2463 = t2461*t2293;
  t2464 = t2404 + t2463;
  t2468 = -0.0115*t2464;
  t2469 = t2261*t2461;
  t2470 = t2469 + t2447;
  t2472 = -0.0265*t2470;
  t2313 = t2260*t2300;
  t2336 = t2323*t2334;
  t2337 = -1.*t2261*t2260;
  t2341 = t2337 + t2338;
  t2342 = -0.0265*t2341;
  t2345 = -1.*t2260*t2293;
  t2346 = t2344 + t2345;
  t2350 = -0.0115*t2346;
  t2351 = t2313 + t2336 + t2342 + t2350;
  t2352 = var1[13]*t2351;
  t2361 = t2141*t2360;
  t2372 = -1.*t2191*t2366;
  t2416 = t2361 + t2372 + t2387 + t2402 + t2408 + t2415;
  t2440 = t2323*t2300;
  t2444 = t2393*t2334;
  t2445 = -0.0115*t2406;
  t2446 = -1.*t2261*t2323;
  t2449 = t2446 + t2447;
  t2450 = -0.0265*t2449;
  t2451 = t2440 + t2444 + t2445 + t2450;
  t2452 = var1[13]*t2451;
  t2453 = -1.*t2191*t2360;
  t2454 = -1.*t2141*t2366;
  t2473 = t2453 + t2454 + t2456 + t2462 + t2468 + t2472;
  t2485 = -0.0695*t2184;
  t2486 = 0.0265*t2222;
  t2488 = t2485 + t2486;
  t2503 = -0.2375*t2261;
  t2505 = 0.0265*t2293;
  t2506 = t2503 + t2505;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=t2352 + (t2387 + t2402 + t2408 + t2415 - 1.*t2191*t2427 - 1.*t2141*t2434)*var1[9] + t2416*var1[12];
  p_output1[5]=t2452 + (-1.*t2141*t2427 + t2191*t2434 + t2456 + t2462 + t2468 + t2472)*var1[9] + t2473*var1[12];
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=t2352 + t2416*var1[9] + (t2361 + t2387 + t2402 + t2408 + t2415 + t2191*t2488)*var1[12];
  p_output1[11]=t2452 + t2473*var1[9] + (t2453 + t2456 + t2462 + t2468 + t2472 + t2141*t2488)*var1[12];
  p_output1[12]=t2351*var1[9] + t2351*var1[12] + (t2313 - 0.0115*(t2345 - 1.*t2261*t2461) - 0.0265*(t2337 + t2463) + t2461*t2506)*var1[13];
  p_output1[13]=t2451*var1[9] + t2451*var1[12] + (-0.0115*t2341 + t2440 - 0.0265*(t2260*t2293 + t2446) + t2260*t2506)*var1[13];
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

#include "dJ_leftToe.hh"

namespace SymFunction
{

void dJ_leftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
