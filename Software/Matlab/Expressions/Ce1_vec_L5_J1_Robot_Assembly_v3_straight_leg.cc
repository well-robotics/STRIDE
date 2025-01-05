/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:27 GMT-05:00
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
  double t473;
  double t430;
  double t436;
  double t486;
  double t498;
  double t466;
  double t491;
  double t492;
  double t429;
  double t505;
  double t507;
  double t512;
  double t516;
  double t526;
  double t527;
  double t529;
  double t495;
  double t517;
  double t537;
  double t538;
  double t539;
  double t519;
  double t545;
  double t546;
  double t547;
  double t548;
  double t549;
  double t550;
  double t551;
  double t552;
  double t531;
  double t533;
  double t520;
  double t521;
  double t557;
  double t570;
  double t571;
  double t572;
  double t573;
  double t574;
  double t575;
  double t576;
  double t585;
  double t586;
  double t581;
  double t582;
  double t587;
  double t588;
  double t589;
  double t590;
  double t592;
  double t593;
  double t594;
  double t595;
  double t599;
  double t600;
  double t601;
  double t602;
  double t603;
  double t583;
  double t584;
  double t591;
  double t596;
  double t597;
  double t607;
  double t608;
  double t609;
  double t598;
  double t604;
  double t605;
  double t618;
  double t619;
  double t620;
  double t621;
  double t622;
  double t615;
  double t616;
  double t613;
  double t614;
  double t624;
  double t636;
  double t637;
  double t639;
  double t640;
  double t641;
  double t627;
  double t629;
  double t645;
  double t646;
  double t647;
  double t638;
  double t642;
  double t643;
  double t631;
  double t659;
  double t660;
  double t661;
  double t662;
  double t663;
  double t665;
  double t666;
  double t667;
  double t668;
  double t669;
  double t673;
  double t674;
  double t697;
  double t698;
  double t699;
  t473 = Cos(var1[5]);
  t430 = Cos(var1[6]);
  t436 = Sin(var1[5]);
  t486 = Sin(var1[6]);
  t498 = Cos(var1[2]);
  t466 = -1.*t430*t436;
  t491 = -1.*t473*t486;
  t492 = t466 + t491;
  t429 = Sin(var1[2]);
  t505 = t473*t430;
  t507 = -1.*t436*t486;
  t512 = t505 + t507;
  t516 = t498*t512;
  t526 = t430*t436;
  t527 = t473*t486;
  t529 = t526 + t527;
  t495 = t429*t492;
  t517 = t495 + t516;
  t537 = t498*t529;
  t538 = t429*t512;
  t539 = t537 + t538;
  t519 = t498*t492;
  t545 = 0.39928*t517*t539;
  t546 = -1.*t473*t430;
  t547 = t436*t486;
  t548 = t546 + t547;
  t549 = t429*t548;
  t550 = t519 + t549;
  t551 = 0.39928*t517*t550;
  t552 = t545 + t551;
  t531 = -1.*t429*t529;
  t533 = t531 + t516;
  t520 = -1.*t429*t512;
  t521 = t519 + t520;
  t557 = -1.*t429*t492;
  t570 = 0.19964*t517*t533;
  t571 = 0.19964*t521*t539;
  t572 = t498*t548;
  t573 = t557 + t572;
  t574 = 0.19964*t517*t573;
  t575 = 0.19964*t521*t550;
  t576 = t570 + t571 + t574 + t575;
  t585 = -1.*t430;
  t586 = 1. + t585;
  t581 = -1.*t473;
  t582 = 1. + t581;
  t587 = -0.0265*t586;
  t588 = -0.025226*t430;
  t589 = -0.07700600000000002*t486;
  t590 = t587 + t588 + t589;
  t592 = -0.2375*t586;
  t593 = -0.314506*t430;
  t594 = -0.0012740000000000008*t486;
  t595 = t592 + t593 + t594;
  t599 = -0.0695*t582;
  t600 = -0.0265*t436;
  t601 = -1.*t436*t590;
  t602 = t473*t595;
  t603 = t599 + t600 + t601 + t602;
  t583 = -0.0265*t582;
  t584 = 0.0695*t436;
  t591 = t473*t590;
  t596 = t436*t595;
  t597 = t583 + t584 + t591 + t596;
  t607 = -1.*t603*t492;
  t608 = -1.*t597*t512;
  t609 = t607 + t608;
  t598 = t597*t529;
  t604 = t603*t512;
  t605 = t598 + t604;
  t618 = -0.0265*t473;
  t619 = -0.0695*t436;
  t620 = -1.*t473*t590;
  t621 = -1.*t436*t595;
  t622 = t618 + t619 + t620 + t621;
  t615 = 0.0695*t473;
  t616 = t615 + t600 + t601 + t602;
  t613 = 0.19964*t517*t609;
  t614 = t603*t492;
  t624 = t597*t512;
  t636 = -0.07700600000000002*t430;
  t637 = t636 + t594;
  t639 = -0.0012740000000000008*t430;
  t640 = 0.07700600000000002*t486;
  t641 = t639 + t640;
  t627 = 0.19964*t605*t550;
  t629 = -1.*t597*t492;
  t645 = -1.*t436*t637;
  t646 = t473*t641;
  t647 = t645 + t646;
  t638 = t473*t637;
  t642 = t436*t641;
  t643 = t638 + t642;
  t631 = -1.*t603*t548;
  t659 = -0.0265*t430;
  t660 = -1.*t430*t590;
  t661 = 0.0695*t486;
  t662 = t595*t486;
  t663 = t659 + t660 + t661 + t662;
  t665 = 0.0695*t430;
  t666 = t430*t595;
  t667 = 0.0265*t486;
  t668 = t590*t486;
  t669 = t665 + t666 + t667 + t668;
  t673 = 0.19964*t663*t517;
  t674 = 0.19964*t669*t550;
  t697 = -0.0002543413600000002*t517;
  t698 = -0.015373477840000005*t550;
  t699 = t697 + t698;
  p_output1[0]=var2[0]*(-0.5*(0.39928*t517*t521 + 0.39928*t533*t539)*var2[2] - 0.5*t552*var2[5] - 0.5*t552*var2[6]);
  p_output1[1]=var2[0]*(-0.5*(0.19964*Power(t521,2) + 0.19964*Power(t533,2) + 0.19964*(t520 - 1.*t498*t529)*t539 + 0.19964*t517*(-1.*t498*t512 + t557))*var2[2] - 0.5*t576*var2[5] - 0.5*t576*var2[6]);
  p_output1[2]=var2[0]*(-0.5*(0.19964*t521*t605 + 0.19964*t533*t609)*var2[2] - 0.5*(t613 + 0.19964*t517*(t614 + t529*t616 + t512*t622 + t624) + t627 + 0.19964*t539*(-1.*t512*t616 - 1.*t492*t622 + t629 + t631))*var2[5] - 0.5*(t613 + t627 + 0.19964*t539*(t629 + t631 - 1.*t512*t643 - 1.*t492*t647) + 0.19964*t517*(t614 + t624 + t529*t643 + t512*t647))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(-0.5*(0.19964*t533*t663 + 0.19964*t521*t669)*var2[2] - 0.5*(t673 + t674)*var2[5] - 0.5*(0.19964*t517*(0.0265*t430 - 0.0695*t486 + t430*t590 - 1.*t486*t595 + t486*t637 + t430*t641) + 0.19964*t539*(-1.*t430*t637 + t486*t641 + t665 + t666 + t667 + t668) + t673 + t674)*var2[6]);
  p_output1[6]=var2[0]*(-0.5*(-0.015373477840000005*t521 - 0.0002543413600000002*t533)*var2[2] - 0.5*t699*var2[5] - 0.5*t699*var2[6]);
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

#include "Ce1_vec_L5_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L5_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
