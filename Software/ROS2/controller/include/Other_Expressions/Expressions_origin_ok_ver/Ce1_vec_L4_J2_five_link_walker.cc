/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:41 GMT-05:00
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
  double t239;
  double t230;
  double t237;
  double t242;
  double t226;
  double t238;
  double t243;
  double t244;
  double t249;
  double t252;
  double t253;
  double t254;
  double t259;
  double t264;
  double t265;
  double t267;
  double t274;
  double t258;
  double t261;
  double t268;
  double t269;
  double t273;
  double t275;
  double t280;
  double t281;
  double t282;
  double t245;
  double t288;
  double t289;
  double t291;
  double t286;
  double t287;
  double t292;
  double t293;
  double t294;
  double t295;
  double t296;
  double t298;
  double t301;
  double t255;
  double t257;
  double t278;
  double t279;
  double t316;
  double t317;
  double t318;
  double t326;
  double t327;
  double t328;
  double t330;
  double t331;
  double t333;
  double t335;
  double t336;
  double t337;
  double t342;
  double t343;
  double t345;
  double t329;
  double t338;
  double t339;
  double t352;
  double t353;
  double t354;
  double t340;
  double t347;
  double t348;
  double t361;
  double t363;
  double t364;
  double t358;
  double t359;
  double t368;
  double t382;
  double t383;
  double t384;
  double t386;
  double t387;
  double t371;
  double t375;
  double t391;
  double t392;
  double t393;
  double t385;
  double t388;
  double t389;
  double t377;
  double t410;
  double t411;
  double t412;
  double t406;
  double t407;
  double t408;
  double t416;
  double t417;
  double t438;
  double t439;
  double t440;
  t239 = Cos(var1[3]);
  t230 = Cos(var1[4]);
  t237 = Sin(var1[3]);
  t242 = Sin(var1[4]);
  t226 = Sin(var1[2]);
  t238 = -1.*t230*t237;
  t243 = -1.*t239*t242;
  t244 = t238 + t243;
  t249 = Cos(var1[2]);
  t252 = t239*t230;
  t253 = -1.*t237*t242;
  t254 = t252 + t253;
  t259 = t249*t254;
  t264 = t230*t237;
  t265 = t239*t242;
  t267 = t264 + t265;
  t274 = -1.*t226*t254;
  t258 = t226*t244;
  t261 = t258 + t259;
  t268 = -1.*t226*t267;
  t269 = t268 + t259;
  t273 = t249*t244;
  t275 = t273 + t274;
  t280 = t249*t267;
  t281 = t226*t254;
  t282 = t280 + t281;
  t245 = -1.*t226*t244;
  t288 = -1.*t239*t230;
  t289 = t237*t242;
  t291 = t288 + t289;
  t286 = 0.14994*t261*t269;
  t287 = 0.14994*t275*t282;
  t292 = t249*t291;
  t293 = t245 + t292;
  t294 = 0.14994*t261*t293;
  t295 = t226*t291;
  t296 = t273 + t295;
  t298 = 0.14994*t275*t296;
  t301 = t286 + t287 + t294 + t298;
  t255 = -1.*t249*t254;
  t257 = t245 + t255;
  t278 = -1.*t249*t267;
  t279 = t278 + t274;
  t316 = 0.29988*t269*t275;
  t317 = 0.29988*t275*t293;
  t318 = t316 + t317;
  t326 = -0.022663*t230;
  t327 = -0.007370999999999989*t242;
  t328 = t326 + t327;
  t330 = -1.*t230;
  t331 = 1. + t330;
  t333 = -0.16*t331;
  t335 = -0.167371*t230;
  t336 = 0.022663*t242;
  t337 = t333 + t335 + t336;
  t342 = -1.*t237*t328;
  t343 = t239*t337;
  t345 = t342 + t343;
  t329 = t239*t328;
  t338 = t237*t337;
  t339 = t329 + t338;
  t352 = -1.*t345*t244;
  t353 = -1.*t339*t254;
  t354 = t352 + t353;
  t340 = t339*t267;
  t347 = t345*t254;
  t348 = t340 + t347;
  t361 = -1.*t239*t328;
  t363 = -1.*t237*t337;
  t364 = t361 + t363;
  t358 = 0.14994*t275*t354;
  t359 = t345*t244;
  t368 = t339*t254;
  t382 = 0.022663*t230;
  t383 = 0.007370999999999989*t242;
  t384 = t382 + t383;
  t386 = -0.007370999999999989*t230;
  t387 = t386 + t336;
  t371 = 0.14994*t348*t293;
  t375 = -1.*t339*t244;
  t391 = t239*t384;
  t392 = -1.*t237*t387;
  t393 = t391 + t392;
  t385 = t237*t384;
  t388 = t239*t387;
  t389 = t385 + t388;
  t377 = -1.*t345*t291;
  t410 = -1.*t230*t328;
  t411 = t337*t242;
  t412 = t410 + t411;
  t406 = t230*t337;
  t407 = t328*t242;
  t408 = t406 + t407;
  t416 = 0.14994*t412*t275;
  t417 = 0.14994*t408*t293;
  t438 = 0.0033980902199999994*t275;
  t439 = -0.0011052077399999983*t293;
  t440 = t438 + t439;
  p_output1[0]=var2[1]*(-0.5*(0.14994*t257*t261 + 0.14994*Power(t269,2) + 0.14994*Power(t275,2) + 0.14994*t279*t282)*var2[2] - 0.5*t301*var2[3] - 0.5*t301*var2[4]);
  p_output1[1]=var2[1]*(-0.5*(0.29988*t257*t275 + 0.29988*t269*t279)*var2[2] - 0.5*t318*var2[3] - 0.5*t318*var2[4]);
  p_output1[2]=var2[1]*(-0.5*(0.14994*t257*t348 + 0.14994*t279*t354)*var2[2] - 0.5*(t358 + 0.14994*t275*(t267*t345 + t359 + t254*t364 + t368) + t371 + 0.14994*t269*(-1.*t254*t345 - 1.*t244*t364 + t375 + t377))*var2[3] - 0.5*(t358 + t371 + 0.14994*t269*(t375 + t377 - 1.*t254*t389 - 1.*t244*t393) + 0.14994*t275*(t359 + t368 + t267*t389 + t254*t393))*var2[4]);
  p_output1[3]=var2[1]*(-0.5*(0.14994*t257*t408 + 0.14994*t279*t412)*var2[2] - 0.5*(t416 + t417)*var2[3] - 0.5*(0.14994*t275*(t230*t328 - 1.*t242*t337 + t230*t384 + t242*t387) + 0.14994*t269*(t242*t384 - 1.*t230*t387 + t406 + t407) + t416 + t417)*var2[4]);
  p_output1[4]=var2[1]*(-0.5*(-0.0011052077399999983*t257 + 0.0033980902199999994*t279)*var2[2] - 0.5*t440*var2[3] - 0.5*t440*var2[4]);
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

#include "Ce1_vec_L4_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L4_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
