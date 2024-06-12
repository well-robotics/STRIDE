/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:39 GMT-05:00
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
  double t217;
  double t208;
  double t209;
  double t223;
  double t227;
  double t213;
  double t224;
  double t225;
  double t205;
  double t230;
  double t231;
  double t234;
  double t237;
  double t244;
  double t245;
  double t246;
  double t226;
  double t238;
  double t252;
  double t253;
  double t254;
  double t239;
  double t259;
  double t260;
  double t261;
  double t262;
  double t263;
  double t264;
  double t265;
  double t266;
  double t248;
  double t249;
  double t240;
  double t242;
  double t271;
  double t284;
  double t285;
  double t286;
  double t287;
  double t288;
  double t289;
  double t290;
  double t295;
  double t296;
  double t297;
  double t299;
  double t300;
  double t301;
  double t302;
  double t303;
  double t304;
  double t308;
  double t309;
  double t310;
  double t298;
  double t305;
  double t306;
  double t314;
  double t315;
  double t316;
  double t307;
  double t311;
  double t312;
  double t323;
  double t324;
  double t325;
  double t320;
  double t321;
  double t327;
  double t339;
  double t340;
  double t341;
  double t343;
  double t344;
  double t330;
  double t332;
  double t348;
  double t349;
  double t350;
  double t342;
  double t345;
  double t346;
  double t334;
  double t362;
  double t363;
  double t364;
  double t366;
  double t367;
  double t368;
  double t372;
  double t373;
  double t394;
  double t395;
  double t396;
  t217 = Cos(var1[3]);
  t208 = Cos(var1[4]);
  t209 = Sin(var1[3]);
  t223 = Sin(var1[4]);
  t227 = Cos(var1[2]);
  t213 = -1.*t208*t209;
  t224 = -1.*t217*t223;
  t225 = t213 + t224;
  t205 = Sin(var1[2]);
  t230 = t217*t208;
  t231 = -1.*t209*t223;
  t234 = t230 + t231;
  t237 = t227*t234;
  t244 = t208*t209;
  t245 = t217*t223;
  t246 = t244 + t245;
  t226 = t205*t225;
  t238 = t226 + t237;
  t252 = t227*t246;
  t253 = t205*t234;
  t254 = t252 + t253;
  t239 = t227*t225;
  t259 = 0.29988*t238*t254;
  t260 = -1.*t217*t208;
  t261 = t209*t223;
  t262 = t260 + t261;
  t263 = t205*t262;
  t264 = t239 + t263;
  t265 = 0.29988*t238*t264;
  t266 = t259 + t265;
  t248 = -1.*t205*t246;
  t249 = t248 + t237;
  t240 = -1.*t205*t234;
  t242 = t239 + t240;
  t271 = -1.*t205*t225;
  t284 = 0.14994*t238*t249;
  t285 = 0.14994*t242*t254;
  t286 = t227*t262;
  t287 = t271 + t286;
  t288 = 0.14994*t238*t287;
  t289 = 0.14994*t242*t264;
  t290 = t284 + t285 + t288 + t289;
  t295 = -0.022663*t208;
  t296 = -0.007370999999999989*t223;
  t297 = t295 + t296;
  t299 = -1.*t208;
  t300 = 1. + t299;
  t301 = -0.16*t300;
  t302 = -0.167371*t208;
  t303 = 0.022663*t223;
  t304 = t301 + t302 + t303;
  t308 = -1.*t209*t297;
  t309 = t217*t304;
  t310 = t308 + t309;
  t298 = t217*t297;
  t305 = t209*t304;
  t306 = t298 + t305;
  t314 = -1.*t310*t225;
  t315 = -1.*t306*t234;
  t316 = t314 + t315;
  t307 = t306*t246;
  t311 = t310*t234;
  t312 = t307 + t311;
  t323 = -1.*t217*t297;
  t324 = -1.*t209*t304;
  t325 = t323 + t324;
  t320 = 0.14994*t238*t316;
  t321 = t310*t225;
  t327 = t306*t234;
  t339 = 0.022663*t208;
  t340 = 0.007370999999999989*t223;
  t341 = t339 + t340;
  t343 = -0.007370999999999989*t208;
  t344 = t343 + t303;
  t330 = 0.14994*t312*t264;
  t332 = -1.*t306*t225;
  t348 = t217*t341;
  t349 = -1.*t209*t344;
  t350 = t348 + t349;
  t342 = t209*t341;
  t345 = t217*t344;
  t346 = t342 + t345;
  t334 = -1.*t310*t262;
  t362 = -1.*t208*t297;
  t363 = t304*t223;
  t364 = t362 + t363;
  t366 = t208*t304;
  t367 = t297*t223;
  t368 = t366 + t367;
  t372 = 0.14994*t364*t238;
  t373 = 0.14994*t368*t264;
  t394 = 0.0033980902199999994*t238;
  t395 = -0.0011052077399999983*t264;
  t396 = t394 + t395;
  p_output1[0]=var2[0]*(-0.5*(0.29988*t238*t242 + 0.29988*t249*t254)*var2[2] - 0.5*t266*var2[3] - 0.5*t266*var2[4]);
  p_output1[1]=var2[0]*(-0.5*(0.14994*Power(t242,2) + 0.14994*Power(t249,2) + 0.14994*(t240 - 1.*t227*t246)*t254 + 0.14994*t238*(-1.*t227*t234 + t271))*var2[2] - 0.5*t290*var2[3] - 0.5*t290*var2[4]);
  p_output1[2]=var2[0]*(-0.5*(0.14994*t242*t312 + 0.14994*t249*t316)*var2[2] - 0.5*(t320 + 0.14994*t238*(t246*t310 + t321 + t234*t325 + t327) + t330 + 0.14994*t254*(-1.*t234*t310 - 1.*t225*t325 + t332 + t334))*var2[3] - 0.5*(t320 + t330 + 0.14994*t254*(t332 + t334 - 1.*t234*t346 - 1.*t225*t350) + 0.14994*t238*(t321 + t327 + t246*t346 + t234*t350))*var2[4]);
  p_output1[3]=var2[0]*(-0.5*(0.14994*t249*t364 + 0.14994*t242*t368)*var2[2] - 0.5*(t372 + t373)*var2[3] - 0.5*(0.14994*t238*(t208*t297 - 1.*t223*t304 + t208*t341 + t223*t344) + 0.14994*t254*(t223*t341 - 1.*t208*t344 + t366 + t367) + t372 + t373)*var2[4]);
  p_output1[4]=var2[0]*(-0.5*(-0.0011052077399999983*t242 + 0.0033980902199999994*t249)*var2[2] - 0.5*t396*var2[3] - 0.5*t396*var2[4]);
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

#include "Ce1_vec_L4_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L4_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
