/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:44 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t261;
  double t257;
  double t258;
  double t263;
  double t279;
  double t280;
  double t281;
  double t283;
  double t284;
  double t285;
  double t286;
  double t287;
  double t289;
  double t272;
  double t273;
  double t275;
  double t270;
  double t294;
  double t295;
  double t296;
  double t252;
  double t302;
  double t303;
  double t305;
  double t259;
  double t264;
  double t265;
  double t282;
  double t292;
  double t293;
  double t313;
  double t317;
  double t319;
  double t320;
  double t327;
  double t329;
  double t298;
  double t307;
  double t309;
  double t269;
  double t338;
  double t340;
  double t343;
  double t355;
  double t356;
  double t357;
  double t333;
  double t335;
  double t348;
  double t382;
  double t383;
  double t385;
  double t388;
  double t390;
  double t358;
  double t359;
  double t360;
  double t365;
  double t368;
  double t369;
  double t371;
  double t397;
  double t398;
  double t399;
  double t386;
  double t391;
  double t394;
  double t376;
  double t276;
  double t277;
  double t336;
  double t347;
  double t351;
  double t414;
  double t312;
  double t316;
  double t370;
  double t374;
  double t378;
  double t425;
  double t395;
  double t400;
  double t401;
  double t427;
  double t428;
  double t429;
  double t403;
  double t404;
  double t405;
  double t450;
  double t451;
  double t452;
  double t454;
  double t455;
  double t456;
  t261 = Cos(var1[3]);
  t257 = Cos(var1[4]);
  t258 = Sin(var1[3]);
  t263 = Sin(var1[4]);
  t279 = -0.022663*t257;
  t280 = -0.007370999999999989*t263;
  t281 = t279 + t280;
  t283 = -1.*t257;
  t284 = 1. + t283;
  t285 = -0.16*t284;
  t286 = -0.167371*t257;
  t287 = 0.022663*t263;
  t289 = t285 + t286 + t287;
  t272 = t261*t257;
  t273 = -1.*t258*t263;
  t275 = t272 + t273;
  t270 = Sin(var1[2]);
  t294 = t257*t258;
  t295 = t261*t263;
  t296 = t294 + t295;
  t252 = Cos(var1[2]);
  t302 = -1.*t258*t281;
  t303 = t261*t289;
  t305 = t302 + t303;
  t259 = -1.*t257*t258;
  t264 = -1.*t261*t263;
  t265 = t259 + t264;
  t282 = t261*t281;
  t292 = t258*t289;
  t293 = t282 + t292;
  t313 = t252*t275;
  t317 = -1.*t305*t265;
  t319 = -1.*t293*t275;
  t320 = t317 + t319;
  t327 = t270*t265;
  t329 = t327 + t313;
  t298 = t293*t296;
  t307 = t305*t275;
  t309 = t298 + t307;
  t269 = t252*t265;
  t338 = -1.*t261*t281;
  t340 = -1.*t258*t289;
  t343 = t338 + t340;
  t355 = -1.*t261*t257;
  t356 = t258*t263;
  t357 = t355 + t356;
  t333 = 0.14994*t329*t320;
  t335 = t305*t265;
  t348 = t293*t275;
  t382 = 0.022663*t257;
  t383 = 0.007370999999999989*t263;
  t385 = t382 + t383;
  t388 = -0.007370999999999989*t257;
  t390 = t388 + t287;
  t358 = t270*t357;
  t359 = t269 + t358;
  t360 = 0.14994*t309*t359;
  t365 = t252*t296;
  t368 = t270*t275;
  t369 = t365 + t368;
  t371 = -1.*t293*t265;
  t397 = t261*t385;
  t398 = -1.*t258*t390;
  t399 = t397 + t398;
  t386 = t258*t385;
  t391 = t261*t390;
  t394 = t386 + t391;
  t376 = -1.*t305*t357;
  t276 = -1.*t270*t275;
  t277 = t269 + t276;
  t336 = t305*t296;
  t347 = t343*t275;
  t351 = t335 + t336 + t347 + t348;
  t414 = -1.*t270*t265;
  t312 = -1.*t270*t296;
  t316 = t312 + t313;
  t370 = -1.*t343*t265;
  t374 = -1.*t305*t275;
  t378 = t370 + t371 + t374 + t376;
  t425 = 0.14994*t277*t320;
  t395 = t394*t296;
  t400 = t399*t275;
  t401 = t335 + t395 + t348 + t400;
  t427 = t252*t357;
  t428 = t414 + t427;
  t429 = 0.14994*t309*t428;
  t403 = -1.*t399*t265;
  t404 = -1.*t394*t275;
  t405 = t371 + t403 + t404 + t376;
  t450 = t257*t289;
  t451 = t281*t263;
  t452 = t450 + t451;
  t454 = -1.*t257*t281;
  t455 = t289*t263;
  t456 = t454 + t455;
  p_output1[0]=var2[2]*(-0.5*(0.14994*t277*t309 + 0.14994*t316*t320)*var2[2] - 0.5*(t333 + 0.14994*t329*t351 + t360 + 0.14994*t369*t378)*var2[3] - 0.5*(t333 + t360 + 0.14994*t329*t401 + 0.14994*t369*t405)*var2[4]);
  p_output1[1]=var2[2]*(-0.5*(0.14994*(t276 - 1.*t252*t296)*t320 + 0.14994*t309*(-1.*t252*t275 + t414))*var2[2] - 0.5*(0.14994*t277*t351 + 0.14994*t316*t378 + t425 + t429)*var2[3] - 0.5*(0.14994*t277*t401 + 0.14994*t316*t405 + t425 + t429)*var2[4]);
  p_output1[2]=var2[2]*(-0.5*(0.29988*t309*t351 + 0.29988*t320*t378)*var2[3] - 0.5*(0.29988*t309*t401 + 0.29988*t320*t405)*var2[4]);
  p_output1[3]=var2[2]*(-0.5*(0.14994*t351*t452 + 0.14994*t378*t456)*var2[3] - 0.5*(0.14994*t309*(t257*t281 - 1.*t263*t289 + t257*t385 + t263*t390) + 0.14994*t320*(t263*t385 - 1.*t257*t390 + t450 + t451) + 0.14994*t401*t452 + 0.14994*t405*t456)*var2[4]);
  p_output1[4]=var2[2]*(-0.5*(-0.0011052077399999983*t351 + 0.0033980902199999994*t378)*var2[3] - 0.5*(-0.0011052077399999983*t401 + 0.0033980902199999994*t405)*var2[4]);
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

#include "Ce1_vec_L4_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L4_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
