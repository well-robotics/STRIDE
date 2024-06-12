/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:33 GMT-05:00
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
  double t263;
  double t269;
  double t241;
  double t293;
  double t315;
  double t339;
  double t267;
  double t295;
  double t314;
  double t350;
  double t354;
  double t355;
  double t360;
  double t389;
  double t392;
  double t402;
  double t403;
  double t405;
  double t366;
  double t369;
  double t382;
  double t424;
  double t425;
  double t393;
  double t399;
  double t401;
  double t407;
  double t409;
  double t411;
  double t412;
  double t414;
  double t415;
  double t417;
  double t419;
  double t420;
  double t421;
  double t327;
  double t341;
  double t342;
  double t356;
  double t363;
  double t376;
  double t378;
  double t379;
  double t381;
  double t384;
  double t450;
  double t451;
  double t452;
  double t442;
  double t426;
  double t427;
  double t429;
  double t431;
  double t432;
  double t449;
  double t453;
  double t454;
  double t455;
  double t456;
  double t457;
  double t458;
  double t459;
  double t349;
  double t364;
  double t365;
  double t367;
  double t368;
  double t370;
  double t371;
  double t372;
  double t373;
  double t374;
  double t380;
  double t388;
  double t422;
  double t438;
  double t439;
  double t440;
  double t441;
  double t443;
  double t444;
  double t445;
  double t446;
  double t447;
  double t448;
  double t460;
  double t468;
  double t469;
  double t470;
  double t482;
  double t483;
  double t484;
  t263 = Cos(var1[5]);
  t269 = Sin(var1[2]);
  t241 = Cos(var1[2]);
  t293 = Sin(var1[5]);
  t315 = Cos(var1[6]);
  t339 = Sin(var1[6]);
  t267 = t241*t263;
  t295 = -1.*t269*t293;
  t314 = t267 + t295;
  t350 = -1.*t263*t269;
  t354 = -1.*t241*t293;
  t355 = t350 + t354;
  t360 = -0.0265*t339;
  t389 = -1.*t315;
  t392 = 1. + t389;
  t402 = -1.*t241*t263;
  t403 = t269*t293;
  t405 = t402 + t403;
  t366 = -1.*t355*t339;
  t369 = t315*t355;
  t382 = -0.0265*t293;
  t424 = -1.*t263;
  t425 = 1. + t424;
  t393 = -0.2375*t392;
  t399 = t393 + t360;
  t401 = t355*t399;
  t407 = -0.0265*t392;
  t409 = 0.2375*t339;
  t411 = t407 + t409;
  t412 = t405*t411;
  t414 = t315*t405;
  t415 = t414 + t366;
  t417 = -0.0265*t415;
  t419 = t405*t339;
  t420 = t369 + t419;
  t421 = 0.0225*t420;
  t327 = -0.0265*t315;
  t341 = -0.2375*t339;
  t342 = t327 + t341;
  t356 = 0.2375*t315;
  t363 = t356 + t360;
  t376 = -0.0265*t263;
  t378 = -0.0695*t293;
  t379 = t376 + t378;
  t381 = 0.0695*t263;
  t384 = t381 + t382;
  t450 = t263*t269;
  t451 = t241*t293;
  t452 = t450 + t451;
  t442 = -1.*t405*t339;
  t426 = -0.0695*t425;
  t427 = t426 + t382;
  t429 = -0.0265*t425;
  t431 = 0.0695*t293;
  t432 = t429 + t431;
  t449 = t405*t399;
  t453 = t452*t411;
  t454 = t452*t339;
  t455 = t414 + t454;
  t456 = 0.0225*t455;
  t457 = t315*t452;
  t458 = t457 + t442;
  t459 = -0.0265*t458;
  t349 = t314*t342;
  t364 = t355*t363;
  t365 = -1.*t315*t314;
  t367 = t365 + t366;
  t368 = -0.0265*t367;
  t370 = -1.*t314*t339;
  t371 = t369 + t370;
  t372 = 0.0225*t371;
  t373 = t349 + t364 + t368 + t372;
  t374 = var1[13]*t373;
  t380 = t241*t379;
  t388 = -1.*t269*t384;
  t422 = t380 + t388 + t401 + t412 + t417 + t421;
  t438 = t355*t342;
  t439 = t405*t363;
  t440 = 0.0225*t415;
  t441 = -1.*t315*t355;
  t443 = t441 + t442;
  t444 = -0.0265*t443;
  t445 = t438 + t439 + t440 + t444;
  t446 = var1[13]*t445;
  t447 = -1.*t269*t379;
  t448 = -1.*t241*t384;
  t460 = t447 + t448 + t449 + t453 + t456 + t459;
  t468 = -0.0695*t263;
  t469 = 0.0265*t293;
  t470 = t468 + t469;
  t482 = -0.2375*t315;
  t483 = 0.0265*t339;
  t484 = t482 + t483;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=t374 + (t401 + t412 + t417 + t421 - 1.*t269*t427 - 1.*t241*t432)*var1[9] + t422*var1[12];
  p_output1[5]=t446 + (-1.*t241*t427 + t269*t432 + t449 + t453 + t456 + t459)*var1[9] + t460*var1[12];
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=t374 + t422*var1[9] + (t380 + t401 + t412 + t417 + t421 + t269*t470)*var1[12];
  p_output1[11]=t446 + t460*var1[9] + (t447 + t449 + t453 + t456 + t459 + t241*t470)*var1[12];
  p_output1[12]=t373*var1[9] + t373*var1[12] + (t349 + 0.0225*(t370 - 1.*t315*t452) - 0.0265*(t365 + t454) + t452*t484)*var1[13];
  p_output1[13]=t445*var1[9] + t445*var1[12] + (0.0225*t367 + t438 - 0.0265*(t314*t339 + t441) + t314*t484)*var1[13];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "dJ_leftToe.hh"

namespace SymFunction
{

void dJ_leftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
