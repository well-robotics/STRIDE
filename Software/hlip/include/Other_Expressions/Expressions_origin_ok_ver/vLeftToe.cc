/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:31 GMT-05:00
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
  double t175;
  double t217;
  double t256;
  double t167;
  double t282;
  double t287;
  double t288;
  double t291;
  double t269;
  double t274;
  double t278;
  double t299;
  double t300;
  double t301;
  double t263;
  double t332;
  double t334;
  double t290;
  double t293;
  double t295;
  double t298;
  double t313;
  double t314;
  double t315;
  double t316;
  double t317;
  double t318;
  double t319;
  double t320;
  double t321;
  double t323;
  double t324;
  double t325;
  double t349;
  double t350;
  double t353;
  double t354;
  double t355;
  double t356;
  double t360;
  double t363;
  double t208;
  double t220;
  double t236;
  double t261;
  double t266;
  double t387;
  double t388;
  double t389;
  double t377;
  double t335;
  double t338;
  double t341;
  double t342;
  double t343;
  double t386;
  double t391;
  double t392;
  double t393;
  double t394;
  double t395;
  double t396;
  double t397;
  t175 = Cos(var1[5]);
  t217 = Sin(var1[5]);
  t256 = Cos(var1[2]);
  t167 = Sin(var1[2]);
  t282 = Cos(var1[6]);
  t287 = -1.*t282;
  t288 = 1. + t287;
  t291 = Sin(var1[6]);
  t269 = t256*t175;
  t274 = -1.*t167*t217;
  t278 = t269 + t274;
  t299 = -1.*t175*t167;
  t300 = -1.*t256*t217;
  t301 = t299 + t300;
  t263 = -0.0265*t217;
  t332 = -1.*t175;
  t334 = 1. + t332;
  t290 = -0.2375*t288;
  t293 = -0.0265*t291;
  t295 = t290 + t293;
  t298 = t278*t295;
  t313 = -0.0265*t288;
  t314 = 0.2375*t291;
  t315 = t313 + t314;
  t316 = t301*t315;
  t317 = t282*t278;
  t318 = t301*t291;
  t319 = t317 + t318;
  t320 = 0.0225*t319;
  t321 = t282*t301;
  t323 = -1.*t278*t291;
  t324 = t321 + t323;
  t325 = -0.0265*t324;
  t349 = t175*t167;
  t350 = t256*t217;
  t353 = t349 + t350;
  t354 = -0.0265*t282;
  t355 = -0.2375*t291;
  t356 = t354 + t355;
  t360 = 0.2375*t282;
  t363 = t360 + t293;
  t208 = -0.0265*t175;
  t220 = -0.0695*t217;
  t236 = t208 + t220;
  t261 = 0.0695*t175;
  t266 = t261 + t263;
  t387 = -1.*t256*t175;
  t388 = t167*t217;
  t389 = t387 + t388;
  t377 = -1.*t301*t291;
  t335 = -0.0695*t334;
  t338 = t335 + t263;
  t341 = -0.0265*t334;
  t342 = 0.0695*t217;
  t343 = t341 + t342;
  t386 = t301*t295;
  t391 = t389*t315;
  t392 = t282*t389;
  t393 = t392 + t377;
  t394 = -0.0265*t393;
  t395 = t389*t291;
  t396 = t321 + t395;
  t397 = 0.0225*t396;
  p_output1[0]=var1[7] + (t298 + t316 + t320 + t325 + t256*t338 - 1.*t167*t343)*var1[9] + (t167*t236 + t256*t266 + t298 + t316 + t320 + t325)*var1[12] + (-0.0265*(t323 - 1.*t282*t353) + 0.0225*(t317 - 1.*t291*t353) + t353*t356 + t278*t363)*var1[13];
  p_output1[1]=var1[8] + (-1.*t167*t338 - 1.*t256*t343 + t386 + t391 + t394 + t397)*var1[9] + (t236*t256 - 1.*t167*t266 + t386 + t391 + t394 + t397)*var1[12] + (0.0225*t324 + t278*t356 + t301*t363 - 0.0265*(-1.*t278*t282 + t377))*var1[13];
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
