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
  double t2057;
  double t2115;
  double t1982;
  double t2159;
  double t2211;
  double t2215;
  double t2218;
  double t2231;
  double t2247;
  double t2248;
  double t2250;
  double t2204;
  double t2207;
  double t2208;
  double t2236;
  double t2294;
  double t2295;
  double t2296;
  double t2269;
  double t2279;
  double t2123;
  double t2324;
  double t2325;
  double t2219;
  double t2237;
  double t2256;
  double t2257;
  double t2265;
  double t2337;
  double t2338;
  double t2339;
  double t2327;
  double t2329;
  double t2332;
  double t2333;
  double t2334;
  double t2268;
  double t2350;
  double t2082;
  double t2133;
  double t2164;
  double t2166;
  double t2179;
  double t2345;
  double t2312;
  double t2315;
  double t2288;
  double t2289;
  double t2298;
  double t2300;
  double t2303;
  double t2346;
  double t2347;
  t2057 = Cos(var1[3]);
  t2115 = Sin(var1[3]);
  t1982 = Sin(var1[2]);
  t2159 = Cos(var1[2]);
  t2211 = Cos(var1[4]);
  t2215 = -1.*t2211;
  t2218 = 1. + t2215;
  t2231 = Sin(var1[4]);
  t2247 = -1.*t2159*t2057;
  t2248 = -1.*t1982*t2115;
  t2250 = t2247 + t2248;
  t2204 = t2057*t1982;
  t2207 = -1.*t2159*t2115;
  t2208 = t2204 + t2207;
  t2236 = -0.2375*t2231;
  t2294 = t2159*t2057;
  t2295 = t1982*t2115;
  t2296 = t2294 + t2295;
  t2269 = -1.*t2208*t2231;
  t2279 = t2211*t2208;
  t2123 = -0.0695*t2115;
  t2324 = -1.*t2057;
  t2325 = 1. + t2324;
  t2219 = -0.0265*t2218;
  t2237 = t2219 + t2236;
  t2256 = -0.2375*t2218;
  t2257 = 0.0265*t2231;
  t2265 = t2256 + t2257;
  t2337 = -1.*t2057*t1982;
  t2338 = t2159*t2115;
  t2339 = t2337 + t2338;
  t2327 = -0.0265*t2325;
  t2329 = t2327 + t2123;
  t2332 = -0.0695*t2325;
  t2333 = 0.0265*t2115;
  t2334 = t2332 + t2333;
  t2268 = t2211*t2250;
  t2350 = t2211*t2339;
  t2082 = 0.0265*t2057;
  t2133 = t2082 + t2123;
  t2164 = -0.0695*t2057;
  t2166 = -0.0265*t2115;
  t2179 = t2164 + t2166;
  t2345 = t2211*t2296;
  t2312 = -1.*t2296*t2231;
  t2315 = t2279 + t2312;
  t2288 = 0.0265*t2211;
  t2289 = t2288 + t2236;
  t2298 = -0.2375*t2211;
  t2300 = -0.0265*t2231;
  t2303 = t2298 + t2300;
  t2346 = -1.*t2339*t2231;
  t2347 = t2345 + t2346;
  p_output1[0]=var1[7] + (t2265*t2296 - 1.*t1982*t2329 + t2159*t2334 + t2237*t2339 + 0.0325*t2347 - 0.0265*(t2231*t2296 + t2350))*var1[9] + (t1982*t2133 + t2159*t2179 + t2208*t2237 + t2250*t2265 + 0.0325*(t2268 + t2269) - 0.0265*(t2231*t2250 + t2279))*var1[10] + (t2208*t2289 + 0.0325*(t2269 - 1.*t2211*t2296) + t2296*t2303 - 0.0265*t2315)*var1[11];
  p_output1[1]=var1[8] + (t2237*t2250 - 1.*t2159*t2329 - 1.*t1982*t2334 + t2265*t2339 - 0.0265*(t2268 + t2231*t2339) + 0.0325*(-1.*t2231*t2250 + t2350))*var1[9] + (t2133*t2159 - 1.*t1982*t2179 + t2208*t2265 + t2237*t2296 + 0.0325*t2315 - 0.0265*(t2208*t2231 + t2345))*var1[10] + (t2289*t2296 + t2303*t2339 + 0.0325*(t2312 - 1.*t2211*t2339) - 0.0265*t2347)*var1[11];
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
