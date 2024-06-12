/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:59 GMT-05:00
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
  double t408;
  double t388;
  double t400;
  double t412;
  double t380;
  double t402;
  double t418;
  double t419;
  double t425;
  double t426;
  double t427;
  double t428;
  double t433;
  double t439;
  double t442;
  double t445;
  double t466;
  double t432;
  double t435;
  double t446;
  double t453;
  double t465;
  double t471;
  double t479;
  double t480;
  double t482;
  double t421;
  double t490;
  double t494;
  double t497;
  double t487;
  double t488;
  double t499;
  double t500;
  double t503;
  double t504;
  double t505;
  double t507;
  double t510;
  double t429;
  double t430;
  double t477;
  double t478;
  double t525;
  double t526;
  double t527;
  double t535;
  double t536;
  double t537;
  double t539;
  double t540;
  double t542;
  double t544;
  double t545;
  double t546;
  double t551;
  double t552;
  double t554;
  double t538;
  double t547;
  double t548;
  double t561;
  double t562;
  double t563;
  double t549;
  double t556;
  double t557;
  double t570;
  double t572;
  double t573;
  double t567;
  double t568;
  double t577;
  double t591;
  double t592;
  double t593;
  double t595;
  double t596;
  double t580;
  double t584;
  double t600;
  double t601;
  double t602;
  double t594;
  double t597;
  double t598;
  double t586;
  double t619;
  double t620;
  double t621;
  double t615;
  double t616;
  double t617;
  double t625;
  double t626;
  double t647;
  double t648;
  double t649;
  t408 = Cos(var1[5]);
  t388 = Cos(var1[6]);
  t400 = Sin(var1[5]);
  t412 = Sin(var1[6]);
  t380 = Sin(var1[2]);
  t402 = -1.*t388*t400;
  t418 = -1.*t408*t412;
  t419 = t402 + t418;
  t425 = Cos(var1[2]);
  t426 = t408*t388;
  t427 = -1.*t400*t412;
  t428 = t426 + t427;
  t433 = t425*t428;
  t439 = t388*t400;
  t442 = t408*t412;
  t445 = t439 + t442;
  t466 = -1.*t380*t428;
  t432 = t380*t419;
  t435 = t432 + t433;
  t446 = -1.*t380*t445;
  t453 = t446 + t433;
  t465 = t425*t419;
  t471 = t465 + t466;
  t479 = t425*t445;
  t480 = t380*t428;
  t482 = t479 + t480;
  t421 = -1.*t380*t419;
  t490 = -1.*t408*t388;
  t494 = t400*t412;
  t497 = t490 + t494;
  t487 = 0.14994*t435*t453;
  t488 = 0.14994*t471*t482;
  t499 = t425*t497;
  t500 = t421 + t499;
  t503 = 0.14994*t435*t500;
  t504 = t380*t497;
  t505 = t465 + t504;
  t507 = 0.14994*t471*t505;
  t510 = t487 + t488 + t503 + t507;
  t429 = -1.*t425*t428;
  t430 = t421 + t429;
  t477 = -1.*t425*t445;
  t478 = t477 + t466;
  t525 = 0.29988*t453*t471;
  t526 = 0.29988*t471*t500;
  t527 = t525 + t526;
  t535 = -0.022659*t388;
  t536 = -0.007367999999999986*t412;
  t537 = t535 + t536;
  t539 = -1.*t388;
  t540 = 1. + t539;
  t542 = -0.16*t540;
  t544 = -0.167368*t388;
  t545 = 0.022659*t412;
  t546 = t542 + t544 + t545;
  t551 = -1.*t400*t537;
  t552 = t408*t546;
  t554 = t551 + t552;
  t538 = t408*t537;
  t547 = t400*t546;
  t548 = t538 + t547;
  t561 = -1.*t554*t419;
  t562 = -1.*t548*t428;
  t563 = t561 + t562;
  t549 = t548*t445;
  t556 = t554*t428;
  t557 = t549 + t556;
  t570 = -1.*t408*t537;
  t572 = -1.*t400*t546;
  t573 = t570 + t572;
  t567 = 0.14994*t471*t563;
  t568 = t554*t419;
  t577 = t548*t428;
  t591 = 0.022659*t388;
  t592 = 0.007367999999999986*t412;
  t593 = t591 + t592;
  t595 = -0.007367999999999986*t388;
  t596 = t595 + t545;
  t580 = 0.14994*t557*t500;
  t584 = -1.*t548*t419;
  t600 = t408*t593;
  t601 = -1.*t400*t596;
  t602 = t600 + t601;
  t594 = t400*t593;
  t597 = t408*t596;
  t598 = t594 + t597;
  t586 = -1.*t554*t497;
  t619 = -1.*t388*t537;
  t620 = t546*t412;
  t621 = t619 + t620;
  t615 = t388*t546;
  t616 = t537*t412;
  t617 = t615 + t616;
  t625 = 0.14994*t621*t471;
  t626 = 0.14994*t617*t500;
  t647 = 0.0033974904599999994*t471;
  t648 = -0.0011047579199999977*t500;
  t649 = t647 + t648;
  p_output1[0]=var2[1]*(-0.5*(0.14994*t430*t435 + 0.14994*Power(t453,2) + 0.14994*Power(t471,2) + 0.14994*t478*t482)*var2[2] - 0.5*t510*var2[5] - 0.5*t510*var2[6]);
  p_output1[1]=var2[1]*(-0.5*(0.29988*t430*t471 + 0.29988*t453*t478)*var2[2] - 0.5*t527*var2[5] - 0.5*t527*var2[6]);
  p_output1[2]=var2[1]*(-0.5*(0.14994*t430*t557 + 0.14994*t478*t563)*var2[2] - 0.5*(t567 + 0.14994*t471*(t445*t554 + t568 + t428*t573 + t577) + t580 + 0.14994*t453*(-1.*t428*t554 - 1.*t419*t573 + t584 + t586))*var2[5] - 0.5*(t567 + t580 + 0.14994*t453*(t584 + t586 - 1.*t428*t598 - 1.*t419*t602) + 0.14994*t471*(t568 + t577 + t445*t598 + t428*t602))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(-0.5*(0.14994*t430*t617 + 0.14994*t478*t621)*var2[2] - 0.5*(t625 + t626)*var2[5] - 0.5*(0.14994*t471*(t388*t537 - 1.*t412*t546 + t388*t593 + t412*t596) + 0.14994*t453*(t412*t593 - 1.*t388*t596 + t615 + t616) + t625 + t626)*var2[6]);
  p_output1[6]=var2[1]*(-0.5*(-0.0011047579199999977*t430 + 0.0033974904599999994*t478)*var2[2] - 0.5*t649*var2[5] - 0.5*t649*var2[6]);
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

#include "Ce1_vec_L5_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L5_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
