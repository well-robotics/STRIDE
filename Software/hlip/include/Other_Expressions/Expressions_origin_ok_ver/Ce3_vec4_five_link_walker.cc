/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:38 GMT-05:00
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
  double t2147;
  double t2201;
  double t2015;
  double t2109;
  double t1822;
  double t2101;
  double t2187;
  double t2191;
  double t2192;
  double t2199;
  double t2200;
  double t2231;
  double t2232;
  double t2241;
  double t2260;
  double t2262;
  double t2266;
  double t2267;
  double t2268;
  double t2272;
  double t2279;
  double t2202;
  double t2210;
  double t2213;
  double t2225;
  double t2228;
  double t2229;
  double t2281;
  double t2286;
  double t2287;
  double t2295;
  double t2296;
  double t2280;
  double t2282;
  double t2297;
  double t2034;
  double t2114;
  double t2121;
  double t2140;
  double t2246;
  double t2312;
  double t2313;
  double t2314;
  double t2372;
  double t2373;
  double t2374;
  double t2351;
  double t2352;
  double t2354;
  double t2359;
  double t2360;
  double t2364;
  double t2396;
  double t2400;
  double t2401;
  double t2403;
  double t2404;
  double t2405;
  double t2230;
  double t2258;
  double t2311;
  double t2315;
  double t2316;
  double t2317;
  double t2320;
  double t2321;
  double t2322;
  double t2323;
  double t2324;
  double t2333;
  double t2371;
  double t2375;
  double t2429;
  double t2430;
  double t2387;
  double t2432;
  double t2433;
  double t2389;
  t2147 = Cos(var1[4]);
  t2201 = Sin(var1[4]);
  t2015 = Sin(var1[2]);
  t2109 = Sin(var1[3]);
  t1822 = Cos(var1[3]);
  t2101 = Cos(var1[2]);
  t2187 = -1.*t2147;
  t2191 = 1. + t2187;
  t2192 = 0.5*t2191;
  t2199 = 0.671885*t2147;
  t2200 = t2192 + t2199;
  t2231 = t1822*t2147;
  t2232 = -1.*t2109*t2201;
  t2241 = t2231 + t2232;
  t2260 = t2200*t2147;
  t2262 = Power(t2201,2);
  t2266 = 0.171885*t2262;
  t2267 = t2260 + t2266;
  t2268 = -1.*t2147*t2109;
  t2272 = -1.*t1822*t2201;
  t2279 = t2268 + t2272;
  t2202 = t2200*t2201;
  t2210 = -0.171885*t2147*t2201;
  t2213 = t2202 + t2210;
  t2225 = t2147*t2109;
  t2228 = t1822*t2201;
  t2229 = t2225 + t2228;
  t2281 = -1.*t2015*t2241;
  t2286 = -1.*t2101*t1822;
  t2287 = t2015*t2109;
  t2295 = t2286 + t2287;
  t2296 = -0.51185934*t2295;
  t2280 = t2101*t2279;
  t2282 = t2280 + t2281;
  t2297 = -1.*t2015*t2279;
  t2034 = -1.*t1822*t2015;
  t2114 = -1.*t2101*t2109;
  t2121 = t2034 + t2114;
  t2140 = -0.51185934*t2121;
  t2246 = t2101*t2241;
  t2312 = -1.*t1822*t2147;
  t2313 = t2109*t2201;
  t2314 = t2312 + t2313;
  t2372 = t1822*t2200;
  t2373 = -0.171885*t2109*t2201;
  t2374 = t2372 + t2373;
  t2351 = -1.*t2200*t2109;
  t2352 = -0.171885*t1822*t2201;
  t2354 = t2351 + t2352;
  t2359 = t2200*t2109;
  t2360 = 0.171885*t1822*t2201;
  t2364 = t2359 + t2360;
  t2396 = -1.*t2200*t2201;
  t2400 = 0.171885*t2147*t2201;
  t2401 = t2396 + t2400;
  t2403 = Power(t2147,2);
  t2404 = -0.171885*t2403;
  t2405 = t2260 + t2404;
  t2230 = -1.*t2015*t2229;
  t2258 = t2230 + t2246;
  t2311 = -0.85216*t2213*t2282;
  t2315 = t2101*t2314;
  t2316 = t2297 + t2315;
  t2317 = -0.85216*t2267*t2316;
  t2320 = t2015*t2279;
  t2321 = t2320 + t2246;
  t2322 = -0.85216*t2213*t2321;
  t2323 = t2015*t2314;
  t2324 = t2280 + t2323;
  t2333 = -0.85216*t2267*t2324;
  t2371 = t2364*t2241;
  t2375 = t2279*t2374;
  t2429 = -0.171885*t2147*t2109;
  t2430 = t2429 + t2352;
  t2387 = -1.*t2279*t2364;
  t2432 = 0.171885*t1822*t2147;
  t2433 = t2432 + t2373;
  t2389 = -1.*t2374*t2314;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(t2140 - 0.85216*t2213*t2258 - 0.85216*t2267*t2282)*var2[0] - 0.5*(-0.85216*t2213*(-1.*t2101*t2229 + t2281) + t2296 - 0.85216*t2267*(-1.*t2101*t2241 + t2297))*var2[1])*var2[3];
  p_output1[3]=(-0.5*(t2140 + t2322 + t2333)*var2[0] - 0.5*(t2296 + t2311 + t2317)*var2[1] - 0.5*(-0.85216*t2267*(t2241*t2354 + t2371 + t2229*t2374 + t2375) - 0.85216*t2213*(-1.*t2279*t2354 - 1.*t2241*t2374 + t2387 + t2389))*var2[2])*var2[3];
  p_output1[4]=var2[3]*(-0.5*(t2322 + t2333 - 0.85216*t2321*t2401 - 0.85216*(t2101*t2229 + t2015*t2241)*t2405)*var2[0] - 0.5*(t2311 + t2317 - 0.85216*t2282*t2401 - 0.85216*t2258*t2405)*var2[1] - 0.5*(-0.85216*(t2229*t2364 + t2241*t2374)*t2401 - 0.85216*(-1.*t2241*t2364 - 1.*t2279*t2374)*t2405 - 0.85216*t2267*(t2371 + t2375 + t2241*t2430 + t2229*t2433) - 0.85216*t2213*(t2387 + t2389 - 1.*t2279*t2430 - 1.*t2241*t2433))*var2[2] - 0.5*(-1.70432*t2267*t2401 - 1.70432*t2213*t2405)*var2[3] + 0.0732367608*t2401*var2[4]);
  p_output1[5]=0;
  p_output1[6]=0;
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

#include "Ce3_vec4_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
