/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:43 GMT-05:00
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
  double t2114;
  double t2176;
  double t2092;
  double t2218;
  double t2268;
  double t2272;
  double t2275;
  double t2288;
  double t2304;
  double t2308;
  double t2309;
  double t2261;
  double t2264;
  double t2265;
  double t2293;
  double t2351;
  double t2352;
  double t2353;
  double t2326;
  double t2336;
  double t2184;
  double t2381;
  double t2382;
  double t2276;
  double t2295;
  double t2313;
  double t2314;
  double t2322;
  double t2394;
  double t2395;
  double t2396;
  double t2384;
  double t2386;
  double t2389;
  double t2390;
  double t2391;
  double t2325;
  double t2407;
  double t2141;
  double t2190;
  double t2222;
  double t2223;
  double t2254;
  double t2402;
  double t2369;
  double t2372;
  double t2345;
  double t2346;
  double t2355;
  double t2357;
  double t2360;
  double t2403;
  double t2404;
  t2114 = Cos(var1[3]);
  t2176 = Sin(var1[3]);
  t2092 = Sin(var1[2]);
  t2218 = Cos(var1[2]);
  t2268 = Cos(var1[4]);
  t2272 = -1.*t2268;
  t2275 = 1. + t2272;
  t2288 = Sin(var1[4]);
  t2304 = -1.*t2218*t2114;
  t2308 = -1.*t2092*t2176;
  t2309 = t2304 + t2308;
  t2261 = t2114*t2092;
  t2264 = -1.*t2218*t2176;
  t2265 = t2261 + t2264;
  t2293 = -0.2375*t2288;
  t2351 = t2218*t2114;
  t2352 = t2092*t2176;
  t2353 = t2351 + t2352;
  t2326 = -1.*t2265*t2288;
  t2336 = t2268*t2265;
  t2184 = -0.0695*t2176;
  t2381 = -1.*t2114;
  t2382 = 1. + t2381;
  t2276 = -0.0265*t2275;
  t2295 = t2276 + t2293;
  t2313 = -0.2375*t2275;
  t2314 = 0.0265*t2288;
  t2322 = t2313 + t2314;
  t2394 = -1.*t2114*t2092;
  t2395 = t2218*t2176;
  t2396 = t2394 + t2395;
  t2384 = -0.0265*t2382;
  t2386 = t2384 + t2184;
  t2389 = -0.0695*t2382;
  t2390 = 0.0265*t2176;
  t2391 = t2389 + t2390;
  t2325 = t2268*t2309;
  t2407 = t2268*t2396;
  t2141 = 0.0265*t2114;
  t2190 = t2141 + t2184;
  t2222 = -0.0695*t2114;
  t2223 = -0.0265*t2176;
  t2254 = t2222 + t2223;
  t2402 = t2268*t2353;
  t2369 = -1.*t2353*t2288;
  t2372 = t2336 + t2369;
  t2345 = 0.0265*t2268;
  t2346 = t2345 + t2293;
  t2355 = -0.2375*t2268;
  t2357 = -0.0265*t2288;
  t2360 = t2355 + t2357;
  t2403 = -1.*t2396*t2288;
  t2404 = t2402 + t2403;
  p_output1[0]=var1[7] + (t2322*t2353 - 1.*t2092*t2386 + t2218*t2391 + t2295*t2396 - 0.0115*t2404 - 0.0265*(t2288*t2353 + t2407))*var1[9] + (t2092*t2190 + t2218*t2254 + t2265*t2295 + t2309*t2322 - 0.0115*(t2325 + t2326) - 0.0265*(t2288*t2309 + t2336))*var1[10] + (t2265*t2346 - 0.0115*(t2326 - 1.*t2268*t2353) + t2353*t2360 - 0.0265*t2372)*var1[11];
  p_output1[1]=var1[8] + (t2295*t2309 - 1.*t2218*t2386 - 1.*t2092*t2391 + t2322*t2396 - 0.0265*(t2325 + t2288*t2396) - 0.0115*(-1.*t2288*t2309 + t2407))*var1[9] + (t2190*t2218 - 1.*t2092*t2254 + t2265*t2322 + t2295*t2353 - 0.0115*t2372 - 0.0265*(t2265*t2288 + t2402))*var1[10] + (t2346*t2353 + t2360*t2396 - 0.0115*(t2369 - 1.*t2268*t2396) - 0.0265*t2404)*var1[11];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "vRightToe.hh"

namespace SymFunction
{

void vRightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
