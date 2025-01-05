/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:41 GMT-05:00
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
  double t2036;
  double t2097;
  double t2161;
  double t2027;
  double t2203;
  double t2208;
  double t2215;
  double t2219;
  double t2191;
  double t2194;
  double t2199;
  double t2231;
  double t2232;
  double t2233;
  double t2184;
  double t2280;
  double t2287;
  double t2218;
  double t2222;
  double t2223;
  double t2226;
  double t2254;
  double t2260;
  double t2261;
  double t2262;
  double t2263;
  double t2264;
  double t2265;
  double t2266;
  double t2268;
  double t2270;
  double t2272;
  double t2274;
  double t2313;
  double t2314;
  double t2315;
  double t2322;
  double t2323;
  double t2325;
  double t2329;
  double t2334;
  double t2092;
  double t2114;
  double t2119;
  double t2176;
  double t2186;
  double t2371;
  double t2372;
  double t2373;
  double t2356;
  double t2288;
  double t2292;
  double t2295;
  double t2300;
  double t2303;
  double t2370;
  double t2374;
  double t2376;
  double t2377;
  double t2380;
  double t2381;
  double t2382;
  double t2383;
  t2036 = Cos(var1[5]);
  t2097 = Sin(var1[5]);
  t2161 = Cos(var1[2]);
  t2027 = Sin(var1[2]);
  t2203 = Cos(var1[6]);
  t2208 = -1.*t2203;
  t2215 = 1. + t2208;
  t2219 = Sin(var1[6]);
  t2191 = t2161*t2036;
  t2194 = -1.*t2027*t2097;
  t2199 = t2191 + t2194;
  t2231 = -1.*t2036*t2027;
  t2232 = -1.*t2161*t2097;
  t2233 = t2231 + t2232;
  t2184 = -0.0265*t2097;
  t2280 = -1.*t2036;
  t2287 = 1. + t2280;
  t2218 = -0.2375*t2215;
  t2222 = -0.0265*t2219;
  t2223 = t2218 + t2222;
  t2226 = t2199*t2223;
  t2254 = -0.0265*t2215;
  t2260 = 0.2375*t2219;
  t2261 = t2254 + t2260;
  t2262 = t2233*t2261;
  t2263 = t2203*t2199;
  t2264 = t2233*t2219;
  t2265 = t2263 + t2264;
  t2266 = -0.0115*t2265;
  t2268 = t2203*t2233;
  t2270 = -1.*t2199*t2219;
  t2272 = t2268 + t2270;
  t2274 = -0.0265*t2272;
  t2313 = t2036*t2027;
  t2314 = t2161*t2097;
  t2315 = t2313 + t2314;
  t2322 = -0.0265*t2203;
  t2323 = -0.2375*t2219;
  t2325 = t2322 + t2323;
  t2329 = 0.2375*t2203;
  t2334 = t2329 + t2222;
  t2092 = -0.0265*t2036;
  t2114 = -0.0695*t2097;
  t2119 = t2092 + t2114;
  t2176 = 0.0695*t2036;
  t2186 = t2176 + t2184;
  t2371 = -1.*t2161*t2036;
  t2372 = t2027*t2097;
  t2373 = t2371 + t2372;
  t2356 = -1.*t2233*t2219;
  t2288 = -0.0695*t2287;
  t2292 = t2288 + t2184;
  t2295 = -0.0265*t2287;
  t2300 = 0.0695*t2097;
  t2303 = t2295 + t2300;
  t2370 = t2233*t2223;
  t2374 = t2373*t2261;
  t2376 = t2203*t2373;
  t2377 = t2376 + t2356;
  t2380 = -0.0265*t2377;
  t2381 = t2373*t2219;
  t2382 = t2268 + t2381;
  t2383 = -0.0115*t2382;
  p_output1[0]=var1[7] + (t2226 + t2262 + t2266 + t2274 + t2161*t2292 - 1.*t2027*t2303)*var1[9] + (t2027*t2119 + t2161*t2186 + t2226 + t2262 + t2266 + t2274)*var1[12] + (-0.0265*(t2270 - 1.*t2203*t2315) - 0.0115*(t2263 - 1.*t2219*t2315) + t2315*t2325 + t2199*t2334)*var1[13];
  p_output1[1]=var1[8] + (-1.*t2027*t2292 - 1.*t2161*t2303 + t2370 + t2374 + t2380 + t2383)*var1[9] + (t2119*t2161 - 1.*t2027*t2186 + t2370 + t2374 + t2380 + t2383)*var1[12] + (-0.0115*t2272 + t2199*t2325 + t2233*t2334 - 0.0265*(-1.*t2199*t2203 + t2356))*var1[13];
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

#include "vLeftToe.hh"

namespace SymFunction
{

void vLeftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
