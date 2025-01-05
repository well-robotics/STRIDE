/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:52 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t160;
  double t146;
  double t149;
  double t168;
  double t182;
  double t158;
  double t169;
  double t175;
  double t125;
  double t191;
  double t192;
  double t197;
  double t199;
  double t212;
  double t214;
  double t218;
  double t176;
  double t203;
  double t223;
  double t224;
  double t225;
  double t231;
  double t232;
  double t234;
  double t235;
  double t237;
  double t238;
  double t239;
  double t241;
  double t242;
  double t221;
  double t222;
  double t204;
  double t208;
  double t209;
  double t261;
  double t262;
  double t263;
  double t264;
  double t265;
  double t266;
  double t267;
  double t276;
  double t277;
  double t272;
  double t273;
  double t278;
  double t279;
  double t280;
  double t281;
  double t283;
  double t284;
  double t285;
  double t286;
  double t290;
  double t291;
  double t292;
  double t293;
  double t294;
  double t274;
  double t275;
  double t282;
  double t287;
  double t288;
  double t298;
  double t299;
  double t300;
  double t289;
  double t295;
  double t296;
  double t308;
  double t309;
  double t310;
  double t311;
  double t312;
  double t305;
  double t306;
  double t317;
  double t304;
  double t314;
  double t327;
  double t328;
  double t329;
  double t331;
  double t332;
  double t318;
  double t320;
  double t336;
  double t337;
  double t338;
  double t330;
  double t333;
  double t334;
  double t322;
  double t356;
  double t357;
  double t358;
  double t359;
  double t360;
  double t350;
  double t351;
  double t352;
  double t353;
  double t354;
  double t364;
  double t365;
  double t388;
  double t389;
  double t390;
  t160 = Cos(var1[3]);
  t146 = Cos(var1[4]);
  t149 = Sin(var1[3]);
  t168 = Sin(var1[4]);
  t182 = Cos(var1[2]);
  t158 = t146*t149;
  t169 = t160*t168;
  t175 = t158 + t169;
  t125 = Sin(var1[2]);
  t191 = t160*t146;
  t192 = -1.*t149*t168;
  t197 = t191 + t192;
  t199 = t182*t197;
  t212 = -1.*t146*t149;
  t214 = -1.*t160*t168;
  t218 = t212 + t214;
  t176 = t125*t175;
  t203 = t176 + t199;
  t223 = t182*t218;
  t224 = t125*t197;
  t225 = t223 + t224;
  t231 = 0.39928*t203*t225;
  t232 = t125*t218;
  t234 = -1.*t160*t146;
  t235 = t149*t168;
  t237 = t234 + t235;
  t238 = t182*t237;
  t239 = t232 + t238;
  t241 = 0.39928*t225*t239;
  t242 = t231 + t241;
  t221 = -1.*t125*t218;
  t222 = t221 + t199;
  t204 = t182*t175;
  t208 = -1.*t125*t197;
  t209 = t204 + t208;
  t261 = 0.19964*t222*t203;
  t262 = 0.19964*t209*t225;
  t263 = 0.19964*t222*t239;
  t264 = -1.*t125*t237;
  t265 = t223 + t264;
  t266 = 0.19964*t225*t265;
  t267 = t261 + t262 + t263 + t266;
  t276 = -1.*t146;
  t277 = 1. + t276;
  t272 = -1.*t160;
  t273 = 1. + t272;
  t278 = -0.2375*t277;
  t279 = -0.314514*t146;
  t280 = 0.0012709999999999978*t168;
  t281 = t278 + t279 + t280;
  t283 = -0.0265*t277;
  t284 = -0.025229*t146;
  t285 = 0.07701400000000003*t168;
  t286 = t283 + t284 + t285;
  t290 = -0.0265*t273;
  t291 = -0.0695*t149;
  t292 = -1.*t149*t281;
  t293 = t160*t286;
  t294 = t290 + t291 + t292 + t293;
  t274 = -0.0695*t273;
  t275 = 0.0265*t149;
  t282 = t160*t281;
  t287 = t149*t286;
  t288 = t274 + t275 + t282 + t287;
  t298 = t294*t218;
  t299 = t288*t197;
  t300 = t298 + t299;
  t289 = -1.*t288*t175;
  t295 = -1.*t294*t197;
  t296 = t289 + t295;
  t308 = -0.0695*t160;
  t309 = -0.0265*t149;
  t310 = -1.*t160*t281;
  t311 = -1.*t149*t286;
  t312 = t308 + t309 + t310 + t311;
  t305 = 0.0265*t160;
  t306 = t305 + t291 + t292 + t293;
  t317 = 0.19964*t225*t300;
  t304 = -1.*t294*t218;
  t314 = -1.*t288*t197;
  t327 = 0.07701400000000003*t146;
  t328 = -0.0012709999999999978*t168;
  t329 = t327 + t328;
  t331 = 0.0012709999999999978*t146;
  t332 = t331 + t285;
  t318 = 0.19964*t296*t239;
  t320 = t288*t218;
  t336 = t160*t329;
  t337 = -1.*t149*t332;
  t338 = t336 + t337;
  t330 = t149*t329;
  t333 = t160*t332;
  t334 = t330 + t333;
  t322 = t294*t237;
  t356 = -0.0695*t146;
  t357 = -1.*t146*t281;
  t358 = 0.0265*t168;
  t359 = t286*t168;
  t360 = t356 + t357 + t358 + t359;
  t350 = 0.0265*t146;
  t351 = t146*t286;
  t352 = 0.0695*t168;
  t353 = t281*t168;
  t354 = t350 + t351 + t352 + t353;
  t364 = 0.19964*t360*t225;
  t365 = 0.19964*t354*t239;
  t388 = 0.015375074960000006*t225;
  t389 = 0.0002537424399999996*t239;
  t390 = t388 + t389;
  p_output1[0]=var2[0]*(-0.5*(0.39928*t203*t209 + 0.39928*t222*t225)*var2[2] - 0.5*t242*var2[3] - 0.5*t242*var2[4]);
  p_output1[1]=var2[0]*(-0.5*(0.19964*(-1.*t125*t175 - 1.*t182*t197)*t203 + 0.19964*Power(t209,2) + 0.19964*Power(t222,2) + 0.19964*(t208 - 1.*t182*t218)*t225)*var2[2] - 0.5*t267*var2[3] - 0.5*t267*var2[4]);
  p_output1[2]=var2[0]*(-0.5*(0.19964*t222*t296 + 0.19964*t209*t300)*var2[2] - 0.5*(0.19964*t225*(t304 - 1.*t175*t306 - 1.*t197*t312 + t314) + t317 + t318 + 0.19964*t203*(t197*t306 + t218*t312 + t320 + t322))*var2[3] - 0.5*(t317 + t318 + 0.19964*t225*(t304 + t314 - 1.*t175*t334 - 1.*t197*t338) + 0.19964*t203*(t320 + t322 + t197*t334 + t218*t338))*var2[4]);
  p_output1[3]=var2[0]*(-0.5*(0.19964*t222*t354 + 0.19964*t209*t360)*var2[2] - 0.5*(t364 + t365)*var2[3] - 0.5*(0.19964*t225*(0.0695*t146 - 0.0265*t168 + t146*t281 - 1.*t168*t286 + t146*t329 + t168*t332) + 0.19964*t203*(t168*t329 - 1.*t146*t332 + t350 + t351 + t352 + t353) + t364 + t365)*var2[4]);
  p_output1[4]=var2[0]*(-0.5*(0.015375074960000006*t209 + 0.0002537424399999996*t222)*var2[2] - 0.5*t390*var2[3] - 0.5*t390*var2[4]);
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

#include "Ce1_vec_L3_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L3_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
