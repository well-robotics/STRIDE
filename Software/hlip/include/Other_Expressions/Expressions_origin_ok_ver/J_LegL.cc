/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:42 GMT-05:00
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
  double t393;
  double t407;
  double t419;
  double t428;
  double t369;
  double t450;
  double t474;
  double t476;
  double t477;
  double t479;
  double t469;
  double t471;
  double t473;
  double t483;
  double t486;
  double t489;
  double t420;
  double t429;
  double t431;
  double t457;
  double t458;
  double t466;
  double t478;
  double t480;
  double t481;
  double t491;
  double t492;
  double t493;
  double t522;
  double t524;
  double t525;
  double t501;
  double t534;
  double t555;
  double t556;
  double t557;
  double t517;
  double t519;
  double t520;
  double t532;
  double t537;
  double t541;
  double t542;
  double t544;
  double t546;
  double t548;
  double t549;
  double t552;
  double t553;
  double t554;
  double t558;
  double t559;
  double t560;
  double t562;
  double t563;
  double t564;
  double t565;
  double t566;
  double t567;
  double t575;
  double t576;
  double t577;
  double t579;
  double t590;
  double t591;
  double t592;
  double t594;
  double t587;
  double t588;
  double t589;
  double t598;
  double t599;
  double t600;
  double t578;
  double t580;
  double t581;
  double t583;
  double t584;
  double t585;
  double t593;
  double t595;
  double t596;
  double t601;
  double t602;
  double t603;
  double t605;
  double t617;
  double t618;
  double t619;
  double t582;
  double t586;
  double t597;
  double t604;
  double t606;
  double t607;
  double t608;
  double t609;
  double t610;
  double t611;
  double t612;
  double t613;
  double t615;
  double t616;
  double t620;
  double t621;
  double t622;
  double t623;
  double t624;
  double t625;
  double t626;
  double t627;
  double t628;
  double t629;
  double t637;
  double t638;
  double t639;
  double t496;
  double t652;
  double t653;
  double t655;
  double t656;
  double t657;
  double t570;
  double t571;
  double t572;
  double t573;
  double t661;
  double t675;
  double t676;
  double t678;
  double t679;
  double t680;
  double t614;
  double t630;
  double t631;
  double t632;
  double t698;
  double t699;
  double t700;
  double t702;
  double t703;
  double t636;
  double t640;
  double t641;
  double t642;
  double t643;
  double t644;
  double t645;
  double t646;
  double t647;
  double t713;
  double t714;
  double t715;
  double t717;
  double t718;
  t393 = Cos(var1[3]);
  t407 = -1.*t393;
  t419 = 1. + t407;
  t428 = Sin(var1[3]);
  t369 = Cos(var1[2]);
  t450 = Sin(var1[2]);
  t474 = Cos(var1[4]);
  t476 = -1.*t474;
  t477 = 1. + t476;
  t479 = Sin(var1[4]);
  t469 = -1.*t369*t393;
  t471 = -1.*t450*t428;
  t473 = t469 + t471;
  t483 = -1.*t393*t450;
  t486 = t369*t428;
  t489 = t483 + t486;
  t420 = -0.0265*t419;
  t429 = -0.0695*t428;
  t431 = t420 + t429;
  t457 = -0.0695*t419;
  t458 = 0.0265*t428;
  t466 = t457 + t458;
  t478 = -0.0265*t477;
  t480 = -0.2375*t479;
  t481 = t478 + t480;
  t491 = -0.2375*t477;
  t492 = 0.0265*t479;
  t493 = t491 + t492;
  t522 = t369*t393;
  t524 = t450*t428;
  t525 = t522 + t524;
  t501 = t474*t489;
  t534 = t474*t525;
  t555 = t393*t450;
  t556 = -1.*t369*t428;
  t557 = t555 + t556;
  t517 = t450*t431;
  t519 = -1.*t369*t466;
  t520 = -1.*t489*t481;
  t532 = -1.*t525*t493;
  t537 = -1.*t489*t479;
  t541 = t534 + t537;
  t542 = -0.0225*t541;
  t544 = t525*t479;
  t546 = t501 + t544;
  t548 = 0.0265*t546;
  t549 = t517 + t519 + t520 + t532 + t542 + t548;
  t552 = -1.*t369*t431;
  t553 = -1.*t450*t466;
  t554 = -1.*t525*t481;
  t558 = -1.*t557*t493;
  t559 = t557*t479;
  t560 = t534 + t559;
  t562 = 0.0265*t560;
  t563 = t474*t557;
  t564 = -1.*t525*t479;
  t565 = t563 + t564;
  t566 = -0.0225*t565;
  t567 = t552 + t553 + t554 + t558 + t562 + t566;
  t575 = Cos(var1[5]);
  t576 = -1.*t575;
  t577 = 1. + t576;
  t579 = Sin(var1[5]);
  t590 = Cos(var1[6]);
  t591 = -1.*t590;
  t592 = 1. + t591;
  t594 = Sin(var1[6]);
  t587 = t369*t575;
  t588 = -1.*t450*t579;
  t589 = t587 + t588;
  t598 = -1.*t575*t450;
  t599 = -1.*t369*t579;
  t600 = t598 + t599;
  t578 = -0.0695*t577;
  t580 = -0.0265*t579;
  t581 = t578 + t580;
  t583 = -0.0265*t577;
  t584 = 0.0695*t579;
  t585 = t583 + t584;
  t593 = -0.2375*t592;
  t595 = -0.0265*t594;
  t596 = t593 + t595;
  t601 = -0.0265*t592;
  t602 = 0.2375*t594;
  t603 = t601 + t602;
  t605 = t590*t589;
  t617 = t575*t450;
  t618 = t369*t579;
  t619 = t617 + t618;
  t582 = -1.*t369*t581;
  t586 = t450*t585;
  t597 = -1.*t589*t596;
  t604 = -1.*t600*t603;
  t606 = t600*t594;
  t607 = t605 + t606;
  t608 = -0.0225*t607;
  t609 = t590*t600;
  t610 = -1.*t589*t594;
  t611 = t609 + t610;
  t612 = 0.0265*t611;
  t613 = t582 + t586 + t597 + t604 + t608 + t612;
  t615 = -1.*t450*t581;
  t616 = -1.*t369*t585;
  t620 = -1.*t619*t596;
  t621 = -1.*t589*t603;
  t622 = -1.*t619*t594;
  t623 = t605 + t622;
  t624 = 0.0265*t623;
  t625 = t590*t619;
  t626 = t589*t594;
  t627 = t625 + t626;
  t628 = -0.0225*t627;
  t629 = t615 + t616 + t620 + t621 + t624 + t628;
  t637 = -1.*t369*t575;
  t638 = t450*t579;
  t639 = t637 + t638;
  t496 = t474*t473;
  t652 = 0.0265*t393;
  t653 = t652 + t429;
  t655 = -0.0695*t393;
  t656 = -0.0265*t428;
  t657 = t655 + t656;
  t570 = Power(t567,2);
  t571 = Power(t549,2);
  t572 = 0.00085849 + t570 + t571;
  t573 = 1/Sqrt(t572);
  t661 = -1.*t557*t479;
  t675 = 0.0265*t474;
  t676 = t675 + t480;
  t678 = -0.2375*t474;
  t679 = -0.0265*t479;
  t680 = t678 + t679;
  t614 = Power(t613,2);
  t630 = Power(t629,2);
  t631 = 0.00085849 + t614 + t630;
  t632 = 1/Sqrt(t631);
  t698 = -0.0265*t575;
  t699 = -0.0695*t579;
  t700 = t698 + t699;
  t702 = 0.0695*t575;
  t703 = t702 + t580;
  t636 = -1.*t600*t596;
  t640 = -1.*t639*t603;
  t641 = t590*t639;
  t642 = -1.*t600*t594;
  t643 = t641 + t642;
  t644 = 0.0265*t643;
  t645 = t639*t594;
  t646 = t609 + t645;
  t647 = -0.0225*t646;
  t713 = -0.0265*t590;
  t714 = -0.2375*t594;
  t715 = t713 + t714;
  t717 = 0.2375*t590;
  t718 = t717 + t595;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0.5*(2.*(t369*t431 + t450*t466 - 1.*t473*t481 - 1.*t489*t493 + 0.0265*(t479*t489 + t496) - 0.0225*(-1.*t473*t479 + t501))*t549 + 2.*t549*t567)*t573;
  p_output1[5]=0.5*t632*(2.*t613*t629 + 2.*t613*(t450*t581 + t369*t585 + t636 + t640 + t644 + t647));
  p_output1[6]=0.5*t573*(2.*t549*(t554 + t558 + t562 + t566 - 1.*t369*t653 + t450*t657) + 2.*t567*(-1.*t473*t493 - 1.*t481*t557 + 0.0265*(t473*t479 + t563) - 1.*t450*t653 - 1.*t369*t657 - 0.0225*(t496 + t661)));
  p_output1[7]=0;
  p_output1[8]=0.5*t573*(2.*t549*(0.0265*t541 - 0.0225*(-1.*t474*t489 + t564) - 1.*t525*t676 - 1.*t489*t680) + 2.*t567*(0.0265*t565 - 0.0225*(-1.*t474*t525 + t661) - 1.*t557*t676 - 1.*t525*t680));
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=0.5*t632*(2.*t629*(t597 + t604 + t608 + t612 - 1.*t450*t700 - 1.*t369*t703) + 2.*t613*(t636 + t640 + t644 + t647 - 1.*t369*t700 + t450*t703));
  p_output1[12]=0;
  p_output1[13]=0.5*t632*(2.*t629*(0.0265*(t610 - 1.*t590*t619) - 0.0225*t623 - 1.*t619*t715 - 1.*t589*t718) + 2.*t613*(-0.0225*t611 + 0.0265*(-1.*t589*t590 + t642) - 1.*t589*t715 - 1.*t600*t718));
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "J_LegL.hh"

namespace SymFunction
{

void J_LegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
