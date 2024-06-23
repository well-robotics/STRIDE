/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:41 GMT-05:00
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
  double t360;
  double t364;
  double t368;
  double t372;
  double t349;
  double t414;
  double t436;
  double t439;
  double t440;
  double t445;
  double t429;
  double t431;
  double t434;
  double t458;
  double t462;
  double t463;
  double t369;
  double t393;
  double t403;
  double t419;
  double t420;
  double t423;
  double t443;
  double t450;
  double t451;
  double t466;
  double t467;
  double t468;
  double t470;
  double t485;
  double t486;
  double t488;
  double t509;
  double t511;
  double t512;
  double t516;
  double t529;
  double t530;
  double t531;
  double t533;
  double t525;
  double t526;
  double t528;
  double t538;
  double t539;
  double t540;
  double t515;
  double t517;
  double t518;
  double t520;
  double t522;
  double t523;
  double t532;
  double t534;
  double t536;
  double t541;
  double t542;
  double t543;
  double t545;
  double t559;
  double t560;
  double t561;
  t360 = Cos(var1[3]);
  t364 = -1.*t360;
  t368 = 1. + t364;
  t372 = Sin(var1[3]);
  t349 = Cos(var1[2]);
  t414 = Sin(var1[2]);
  t436 = Cos(var1[4]);
  t439 = -1.*t436;
  t440 = 1. + t439;
  t445 = Sin(var1[4]);
  t429 = t349*t360;
  t431 = t414*t372;
  t434 = t429 + t431;
  t458 = t360*t414;
  t462 = -1.*t349*t372;
  t463 = t458 + t462;
  t369 = -0.0265*t368;
  t393 = -0.0695*t372;
  t403 = t369 + t393;
  t419 = -0.0695*t368;
  t420 = 0.0265*t372;
  t423 = t419 + t420;
  t443 = -0.0265*t440;
  t450 = -0.2375*t445;
  t451 = t443 + t450;
  t466 = -0.2375*t440;
  t467 = 0.0265*t445;
  t468 = t466 + t467;
  t470 = t436*t434;
  t485 = -1.*t360*t414;
  t486 = t349*t372;
  t488 = t485 + t486;
  t509 = Cos(var1[5]);
  t511 = -1.*t509;
  t512 = 1. + t511;
  t516 = Sin(var1[5]);
  t529 = Cos(var1[6]);
  t530 = -1.*t529;
  t531 = 1. + t530;
  t533 = Sin(var1[6]);
  t525 = t349*t509;
  t526 = -1.*t414*t516;
  t528 = t525 + t526;
  t538 = -1.*t509*t414;
  t539 = -1.*t349*t516;
  t540 = t538 + t539;
  t515 = -0.0695*t512;
  t517 = -0.0265*t516;
  t518 = t515 + t517;
  t520 = -0.0265*t512;
  t522 = 0.0695*t516;
  t523 = t520 + t522;
  t532 = -0.2375*t531;
  t534 = -0.0265*t533;
  t536 = t532 + t534;
  t541 = -0.0265*t531;
  t542 = 0.2375*t533;
  t543 = t541 + t542;
  t545 = t529*t528;
  t559 = t509*t414;
  t560 = t349*t516;
  t561 = t559 + t560;
  p_output1[0]=Sqrt(0.00085849 + Power(-1.*t349*t403 - 1.*t414*t423 - 1.*t434*t451 - 0.0225*(-1.*t434*t445 + t436*t463) - 1.*t463*t468 + 0.0265*(t445*t463 + t470),2) + Power(t403*t414 - 1.*t349*t423 - 1.*t434*t468 - 1.*t451*t488 + 0.0265*(t434*t445 + t436*t488) - 0.0225*(t470 - 1.*t445*t488),2));
  p_output1[1]=Sqrt(0.00085849 + Power(-1.*t349*t518 + t414*t523 - 1.*t528*t536 + 0.0265*(-1.*t528*t533 + t529*t540) - 1.*t540*t543 - 0.0225*(t533*t540 + t545),2) + Power(-1.*t414*t518 - 1.*t349*t523 - 1.*t528*t543 - 1.*t536*t561 - 0.0225*(t528*t533 + t529*t561) + 0.0265*(t545 - 1.*t533*t561),2));
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
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
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

#include "LegL.hh"

namespace SymFunction
{

void LegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
