/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:45 GMT-05:00
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
  double t2123;
  double t2134;
  double t2082;
  double t2164;
  double t2204;
  double t2236;
  double t2133;
  double t2166;
  double t2203;
  double t2257;
  double t2265;
  double t2266;
  double t2272;
  double t2316;
  double t2319;
  double t2332;
  double t2333;
  double t2336;
  double t2281;
  double t2287;
  double t2308;
  double t2364;
  double t2366;
  double t2320;
  double t2327;
  double t2330;
  double t2338;
  double t2341;
  double t2342;
  double t2345;
  double t2347;
  double t2349;
  double t2351;
  double t2352;
  double t2356;
  double t2358;
  double t2219;
  double t2237;
  double t2243;
  double t2268;
  double t2273;
  double t2298;
  double t2300;
  double t2303;
  double t2307;
  double t2309;
  double t2400;
  double t2401;
  double t2404;
  double t2390;
  double t2367;
  double t2370;
  double t2374;
  double t2375;
  double t2377;
  double t2399;
  double t2405;
  double t2406;
  double t2407;
  double t2411;
  double t2412;
  double t2413;
  double t2415;
  double t2256;
  double t2279;
  double t2280;
  double t2284;
  double t2285;
  double t2288;
  double t2289;
  double t2293;
  double t2294;
  double t2295;
  double t2304;
  double t2315;
  double t2359;
  double t2383;
  double t2387;
  double t2388;
  double t2389;
  double t2392;
  double t2393;
  double t2394;
  double t2395;
  double t2396;
  double t2397;
  double t2416;
  double t2428;
  double t2429;
  double t2431;
  double t2446;
  double t2448;
  double t2449;
  t2123 = Cos(var1[5]);
  t2134 = Sin(var1[2]);
  t2082 = Cos(var1[2]);
  t2164 = Sin(var1[5]);
  t2204 = Cos(var1[6]);
  t2236 = Sin(var1[6]);
  t2133 = t2082*t2123;
  t2166 = -1.*t2134*t2164;
  t2203 = t2133 + t2166;
  t2257 = -1.*t2123*t2134;
  t2265 = -1.*t2082*t2164;
  t2266 = t2257 + t2265;
  t2272 = -0.0265*t2236;
  t2316 = -1.*t2204;
  t2319 = 1. + t2316;
  t2332 = -1.*t2082*t2123;
  t2333 = t2134*t2164;
  t2336 = t2332 + t2333;
  t2281 = -1.*t2266*t2236;
  t2287 = t2204*t2266;
  t2308 = -0.0265*t2164;
  t2364 = -1.*t2123;
  t2366 = 1. + t2364;
  t2320 = -0.2375*t2319;
  t2327 = t2320 + t2272;
  t2330 = t2266*t2327;
  t2338 = -0.0265*t2319;
  t2341 = 0.2375*t2236;
  t2342 = t2338 + t2341;
  t2345 = t2336*t2342;
  t2347 = t2204*t2336;
  t2349 = t2347 + t2281;
  t2351 = -0.0265*t2349;
  t2352 = t2336*t2236;
  t2356 = t2287 + t2352;
  t2358 = 0.0325*t2356;
  t2219 = -0.0265*t2204;
  t2237 = -0.2375*t2236;
  t2243 = t2219 + t2237;
  t2268 = 0.2375*t2204;
  t2273 = t2268 + t2272;
  t2298 = -0.0265*t2123;
  t2300 = -0.0695*t2164;
  t2303 = t2298 + t2300;
  t2307 = 0.0695*t2123;
  t2309 = t2307 + t2308;
  t2400 = t2123*t2134;
  t2401 = t2082*t2164;
  t2404 = t2400 + t2401;
  t2390 = -1.*t2336*t2236;
  t2367 = -0.0695*t2366;
  t2370 = t2367 + t2308;
  t2374 = -0.0265*t2366;
  t2375 = 0.0695*t2164;
  t2377 = t2374 + t2375;
  t2399 = t2336*t2327;
  t2405 = t2404*t2342;
  t2406 = t2404*t2236;
  t2407 = t2347 + t2406;
  t2411 = 0.0325*t2407;
  t2412 = t2204*t2404;
  t2413 = t2412 + t2390;
  t2415 = -0.0265*t2413;
  t2256 = t2203*t2243;
  t2279 = t2266*t2273;
  t2280 = -1.*t2204*t2203;
  t2284 = t2280 + t2281;
  t2285 = -0.0265*t2284;
  t2288 = -1.*t2203*t2236;
  t2289 = t2287 + t2288;
  t2293 = 0.0325*t2289;
  t2294 = t2256 + t2279 + t2285 + t2293;
  t2295 = var1[13]*t2294;
  t2304 = t2082*t2303;
  t2315 = -1.*t2134*t2309;
  t2359 = t2304 + t2315 + t2330 + t2345 + t2351 + t2358;
  t2383 = t2266*t2243;
  t2387 = t2336*t2273;
  t2388 = 0.0325*t2349;
  t2389 = -1.*t2204*t2266;
  t2392 = t2389 + t2390;
  t2393 = -0.0265*t2392;
  t2394 = t2383 + t2387 + t2388 + t2393;
  t2395 = var1[13]*t2394;
  t2396 = -1.*t2134*t2303;
  t2397 = -1.*t2082*t2309;
  t2416 = t2396 + t2397 + t2399 + t2405 + t2411 + t2415;
  t2428 = -0.0695*t2123;
  t2429 = 0.0265*t2164;
  t2431 = t2428 + t2429;
  t2446 = -0.2375*t2204;
  t2448 = 0.0265*t2236;
  t2449 = t2446 + t2448;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=t2295 + (t2330 + t2345 + t2351 + t2358 - 1.*t2134*t2370 - 1.*t2082*t2377)*var1[9] + t2359*var1[12];
  p_output1[5]=t2395 + (-1.*t2082*t2370 + t2134*t2377 + t2399 + t2405 + t2411 + t2415)*var1[9] + t2416*var1[12];
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=t2295 + t2359*var1[9] + (t2304 + t2330 + t2345 + t2351 + t2358 + t2134*t2431)*var1[12];
  p_output1[11]=t2395 + t2416*var1[9] + (t2396 + t2399 + t2405 + t2411 + t2415 + t2082*t2431)*var1[12];
  p_output1[12]=t2294*var1[9] + t2294*var1[12] + (t2256 + 0.0325*(t2288 - 1.*t2204*t2404) - 0.0265*(t2280 + t2406) + t2404*t2449)*var1[13];
  p_output1[13]=t2394*var1[9] + t2394*var1[12] + (0.0325*t2284 + t2383 - 0.0265*(t2203*t2236 + t2389) + t2203*t2449)*var1[13];
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
