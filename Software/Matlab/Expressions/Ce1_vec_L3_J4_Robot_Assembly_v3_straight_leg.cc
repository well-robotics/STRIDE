/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:59 GMT-05:00
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
  double t247;
  double t261;
  double t252;
  double t256;
  double t301;
  double t299;
  double t271;
  double t278;
  double t282;
  double t283;
  double t257;
  double t260;
  double t262;
  double t266;
  double t316;
  double t293;
  double t318;
  double t319;
  double t320;
  double t330;
  double t339;
  double t340;
  double t348;
  double t349;
  double t300;
  double t307;
  double t313;
  double t251;
  double t269;
  double t270;
  double t287;
  double t290;
  double t354;
  double t359;
  double t362;
  double t321;
  double t372;
  double t373;
  double t374;
  double t375;
  double t402;
  double t403;
  double t404;
  double t398;
  double t399;
  double t376;
  double t378;
  double t379;
  double t381;
  double t383;
  double t387;
  double t388;
  double t367;
  double t315;
  double t323;
  double t438;
  double t416;
  double t417;
  double t418;
  double t419;
  double t420;
  double t421;
  double t422;
  double t401;
  double t409;
  double t410;
  double t366;
  double t368;
  double t440;
  double t441;
  double t442;
  double t459;
  double t460;
  double t461;
  double t455;
  double t457;
  double t467;
  double t468;
  double t469;
  double t470;
  double t471;
  double t473;
  double t474;
  double t475;
  double t476;
  double t477;
  double t464;
  double t465;
  double t458;
  double t462;
  double t463;
  double t478;
  double t482;
  double t501;
  double t502;
  double t503;
  double t497;
  double t498;
  double t499;
  double t484;
  t247 = Cos(var1[4]);
  t261 = Sin(var1[4]);
  t252 = -1.*t247;
  t256 = 1. + t252;
  t301 = Cos(var1[3]);
  t299 = Sin(var1[3]);
  t271 = -0.2375*t256;
  t278 = -0.314514*t247;
  t282 = 0.0012709999999999978*t261;
  t283 = t271 + t278 + t282;
  t257 = -0.0265*t256;
  t260 = -0.025229*t247;
  t262 = 0.07701400000000003*t261;
  t266 = t257 + t260 + t262;
  t316 = Cos(var1[2]);
  t293 = Sin(var1[2]);
  t318 = t301*t247;
  t319 = -1.*t299*t261;
  t320 = t318 + t319;
  t330 = -0.0695*t247;
  t339 = -1.*t247*t283;
  t340 = 0.0265*t261;
  t348 = t266*t261;
  t349 = t330 + t339 + t340 + t348;
  t300 = -1.*t247*t299;
  t307 = -1.*t301*t261;
  t313 = t300 + t307;
  t251 = 0.0265*t247;
  t269 = t247*t266;
  t270 = 0.0695*t261;
  t287 = t283*t261;
  t290 = t251 + t269 + t270 + t287;
  t354 = t247*t299;
  t359 = t301*t261;
  t362 = t354 + t359;
  t321 = t316*t320;
  t372 = t316*t313;
  t373 = t293*t320;
  t374 = t372 + t373;
  t375 = 0.19964*t349*t374;
  t402 = 0.07701400000000003*t247;
  t403 = -0.0012709999999999978*t261;
  t404 = t402 + t403;
  t398 = 0.0012709999999999978*t247;
  t399 = t398 + t262;
  t376 = t293*t313;
  t378 = -1.*t301*t247;
  t379 = t299*t261;
  t381 = t378 + t379;
  t383 = t316*t381;
  t387 = t376 + t383;
  t388 = 0.19964*t290*t387;
  t367 = -1.*t293*t320;
  t315 = -1.*t293*t313;
  t323 = t315 + t321;
  t438 = 0.19964*t349*t323;
  t416 = 0.0695*t247;
  t417 = t247*t404;
  t418 = t247*t283;
  t419 = -0.0265*t261;
  t420 = -1.*t266*t261;
  t421 = t399*t261;
  t422 = t416 + t417 + t418 + t419 + t420 + t421;
  t401 = -1.*t247*t399;
  t409 = t404*t261;
  t410 = t251 + t269 + t401 + t270 + t409 + t287;
  t366 = t316*t362;
  t368 = t366 + t367;
  t440 = -1.*t293*t381;
  t441 = t372 + t440;
  t442 = 0.19964*t290*t441;
  t459 = -0.0695*t299;
  t460 = -1.*t299*t283;
  t461 = t301*t266;
  t455 = -1.*t301;
  t457 = 1. + t455;
  t467 = -0.0695*t301;
  t468 = -0.0265*t299;
  t469 = -1.*t301*t283;
  t470 = -1.*t299*t266;
  t471 = t467 + t468 + t469 + t470;
  t473 = -0.0695*t457;
  t474 = 0.0265*t299;
  t475 = t301*t283;
  t476 = t299*t266;
  t477 = t473 + t474 + t475 + t476;
  t464 = 0.0265*t301;
  t465 = t464 + t459 + t460 + t461;
  t458 = -0.0265*t457;
  t462 = t458 + t459 + t460 + t461;
  t463 = -1.*t462*t313;
  t478 = -1.*t477*t320;
  t482 = t477*t313;
  t501 = t301*t404;
  t502 = -1.*t299*t399;
  t503 = t501 + t502;
  t497 = t299*t404;
  t498 = t301*t399;
  t499 = t497 + t498;
  t484 = t462*t381;
  p_output1[0]=var2[3]*(-0.5*(0.19964*t290*t323 + 0.19964*t349*t368)*var2[2] - 0.5*(t375 + t388)*var2[3] - 0.5*(t375 + t388 + 0.19964*(t321 + t293*t362)*t410 + 0.19964*t374*t422)*var2[4]);
  p_output1[1]=var2[3]*(-0.5*(0.19964*t349*(-1.*t316*t320 - 1.*t293*t362) + 0.19964*t290*(-1.*t313*t316 + t367))*var2[2] - 0.5*(t438 + t442)*var2[3] - 0.5*(0.19964*t368*t410 + 0.19964*t323*t422 + t438 + t442)*var2[4]);
  p_output1[2]=var2[3]*(-0.5*(0.19964*t290*(t463 - 1.*t362*t465 - 1.*t320*t471 + t478) + 0.19964*t349*(t320*t465 + t313*t471 + t482 + t484))*var2[3] - 0.5*(0.19964*t410*(t313*t462 + t320*t477) + 0.19964*t422*(-1.*t320*t462 - 1.*t362*t477) + 0.19964*t349*(t482 + t484 + t320*t499 + t313*t503) + 0.19964*t290*(t463 + t478 - 1.*t362*t499 - 1.*t320*t503))*var2[4]);
  p_output1[3]=-0.5*(0.39928*t349*t410 + 0.39928*t290*t422)*var2[3]*var2[4];
  p_output1[4]=-0.5*(0.015375074960000006*t410 + 0.0002537424399999996*t422)*var2[3]*var2[4];
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

#include "Ce1_vec_L3_J4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L3_J4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
