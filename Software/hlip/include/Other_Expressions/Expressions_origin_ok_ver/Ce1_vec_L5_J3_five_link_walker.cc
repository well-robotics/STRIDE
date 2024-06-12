/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:38:01 GMT-05:00
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
  double t435;
  double t430;
  double t432;
  double t438;
  double t478;
  double t479;
  double t480;
  double t484;
  double t485;
  double t486;
  double t487;
  double t488;
  double t494;
  double t460;
  double t465;
  double t471;
  double t455;
  double t503;
  double t504;
  double t505;
  double t426;
  double t511;
  double t512;
  double t514;
  double t433;
  double t439;
  double t442;
  double t482;
  double t499;
  double t500;
  double t522;
  double t526;
  double t528;
  double t529;
  double t536;
  double t538;
  double t507;
  double t516;
  double t518;
  double t453;
  double t547;
  double t549;
  double t552;
  double t564;
  double t565;
  double t566;
  double t542;
  double t544;
  double t557;
  double t591;
  double t592;
  double t594;
  double t597;
  double t599;
  double t567;
  double t568;
  double t569;
  double t574;
  double t577;
  double t578;
  double t580;
  double t606;
  double t607;
  double t608;
  double t595;
  double t600;
  double t603;
  double t585;
  double t473;
  double t475;
  double t545;
  double t556;
  double t560;
  double t623;
  double t521;
  double t525;
  double t579;
  double t583;
  double t587;
  double t634;
  double t604;
  double t609;
  double t610;
  double t636;
  double t637;
  double t638;
  double t612;
  double t613;
  double t614;
  double t659;
  double t660;
  double t661;
  double t663;
  double t664;
  double t665;
  t435 = Cos(var1[5]);
  t430 = Cos(var1[6]);
  t432 = Sin(var1[5]);
  t438 = Sin(var1[6]);
  t478 = -0.022659*t430;
  t479 = -0.007367999999999986*t438;
  t480 = t478 + t479;
  t484 = -1.*t430;
  t485 = 1. + t484;
  t486 = -0.16*t485;
  t487 = -0.167368*t430;
  t488 = 0.022659*t438;
  t494 = t486 + t487 + t488;
  t460 = t435*t430;
  t465 = -1.*t432*t438;
  t471 = t460 + t465;
  t455 = Sin(var1[2]);
  t503 = t430*t432;
  t504 = t435*t438;
  t505 = t503 + t504;
  t426 = Cos(var1[2]);
  t511 = -1.*t432*t480;
  t512 = t435*t494;
  t514 = t511 + t512;
  t433 = -1.*t430*t432;
  t439 = -1.*t435*t438;
  t442 = t433 + t439;
  t482 = t435*t480;
  t499 = t432*t494;
  t500 = t482 + t499;
  t522 = t426*t471;
  t526 = -1.*t514*t442;
  t528 = -1.*t500*t471;
  t529 = t526 + t528;
  t536 = t455*t442;
  t538 = t536 + t522;
  t507 = t500*t505;
  t516 = t514*t471;
  t518 = t507 + t516;
  t453 = t426*t442;
  t547 = -1.*t435*t480;
  t549 = -1.*t432*t494;
  t552 = t547 + t549;
  t564 = -1.*t435*t430;
  t565 = t432*t438;
  t566 = t564 + t565;
  t542 = 0.14994*t538*t529;
  t544 = t514*t442;
  t557 = t500*t471;
  t591 = 0.022659*t430;
  t592 = 0.007367999999999986*t438;
  t594 = t591 + t592;
  t597 = -0.007367999999999986*t430;
  t599 = t597 + t488;
  t567 = t455*t566;
  t568 = t453 + t567;
  t569 = 0.14994*t518*t568;
  t574 = t426*t505;
  t577 = t455*t471;
  t578 = t574 + t577;
  t580 = -1.*t500*t442;
  t606 = t435*t594;
  t607 = -1.*t432*t599;
  t608 = t606 + t607;
  t595 = t432*t594;
  t600 = t435*t599;
  t603 = t595 + t600;
  t585 = -1.*t514*t566;
  t473 = -1.*t455*t471;
  t475 = t453 + t473;
  t545 = t514*t505;
  t556 = t552*t471;
  t560 = t544 + t545 + t556 + t557;
  t623 = -1.*t455*t442;
  t521 = -1.*t455*t505;
  t525 = t521 + t522;
  t579 = -1.*t552*t442;
  t583 = -1.*t514*t471;
  t587 = t579 + t580 + t583 + t585;
  t634 = 0.14994*t475*t529;
  t604 = t603*t505;
  t609 = t608*t471;
  t610 = t544 + t604 + t557 + t609;
  t636 = t426*t566;
  t637 = t623 + t636;
  t638 = 0.14994*t518*t637;
  t612 = -1.*t608*t442;
  t613 = -1.*t603*t471;
  t614 = t580 + t612 + t613 + t585;
  t659 = t430*t494;
  t660 = t480*t438;
  t661 = t659 + t660;
  t663 = -1.*t430*t480;
  t664 = t494*t438;
  t665 = t663 + t664;
  p_output1[0]=var2[2]*(-0.5*(0.14994*t475*t518 + 0.14994*t525*t529)*var2[2] - 0.5*(t542 + 0.14994*t538*t560 + t569 + 0.14994*t578*t587)*var2[5] - 0.5*(t542 + t569 + 0.14994*t538*t610 + 0.14994*t578*t614)*var2[6]);
  p_output1[1]=var2[2]*(-0.5*(0.14994*(t473 - 1.*t426*t505)*t529 + 0.14994*t518*(-1.*t426*t471 + t623))*var2[2] - 0.5*(0.14994*t475*t560 + 0.14994*t525*t587 + t634 + t638)*var2[5] - 0.5*(0.14994*t475*t610 + 0.14994*t525*t614 + t634 + t638)*var2[6]);
  p_output1[2]=var2[2]*(-0.5*(0.29988*t518*t560 + 0.29988*t529*t587)*var2[5] - 0.5*(0.29988*t518*t610 + 0.29988*t529*t614)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[2]*(-0.5*(0.14994*t560*t661 + 0.14994*t587*t665)*var2[5] - 0.5*(0.14994*t518*(t430*t480 - 1.*t438*t494 + t430*t594 + t438*t599) + 0.14994*t529*(t438*t594 - 1.*t430*t599 + t659 + t660) + 0.14994*t610*t661 + 0.14994*t614*t665)*var2[6]);
  p_output1[6]=var2[2]*(-0.5*(-0.0011047579199999977*t560 + 0.0033974904599999994*t587)*var2[5] - 0.5*(-0.0011047579199999977*t610 + 0.0033974904599999994*t614)*var2[6]);
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

#include "Ce1_vec_L5_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L5_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
