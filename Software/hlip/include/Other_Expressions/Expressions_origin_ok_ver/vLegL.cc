/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:43 GMT-05:00
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
  double t429;
  double t457;
  double t420;
  double t480;
  double t498;
  double t500;
  double t502;
  double t508;
  double t524;
  double t534;
  double t541;
  double t494;
  double t496;
  double t497;
  double t458;
  double t571;
  double t572;
  double t506;
  double t515;
  double t517;
  double t542;
  double t544;
  double t546;
  double t587;
  double t593;
  double t595;
  double t563;
  double t443;
  double t467;
  double t482;
  double t486;
  double t491;
  double t601;
  double t602;
  double t605;
  double t606;
  double t607;
  double t609;
  double t611;
  double t614;
  double t617;
  double t573;
  double t574;
  double t580;
  double t583;
  double t584;
  double t630;
  double t631;
  double t633;
  double t578;
  double t586;
  double t618;
  double t627;
  double t628;
  double t634;
  double t635;
  double t638;
  double t641;
  double t643;
  double t645;
  double t646;
  double t648;
  double t649;
  double t650;
  double t550;
  double t654;
  double t655;
  double t656;
  double t658;
  double t555;
  double t677;
  double t678;
  double t681;
  double t682;
  double t683;
  double t704;
  double t706;
  double t709;
  double t711;
  double t705;
  double t707;
  double t708;
  double t716;
  double t717;
  double t719;
  double t733;
  double t734;
  double t721;
  double t743;
  double t744;
  double t728;
  double t729;
  double t730;
  double t710;
  double t712;
  double t713;
  double t720;
  double t722;
  double t752;
  double t759;
  double t760;
  double t761;
  double t735;
  double t736;
  double t737;
  double t739;
  double t740;
  double t741;
  double t745;
  double t746;
  double t748;
  double t749;
  double t750;
  double t764;
  double t765;
  double t738;
  double t742;
  double t747;
  double t751;
  double t753;
  double t754;
  double t755;
  double t756;
  double t757;
  double t771;
  double t772;
  double t773;
  double t774;
  double t775;
  double t776;
  double t777;
  double t778;
  double t779;
  double t780;
  double t783;
  double t784;
  double t785;
  double t786;
  double t788;
  double t789;
  double t790;
  double t792;
  double t793;
  double t800;
  double t801;
  double t802;
  double t725;
  double t799;
  double t803;
  double t804;
  double t805;
  double t806;
  double t807;
  double t808;
  double t809;
  t429 = Cos(var1[3]);
  t457 = Sin(var1[3]);
  t420 = Sin(var1[2]);
  t480 = Cos(var1[2]);
  t498 = Cos(var1[4]);
  t500 = -1.*t498;
  t502 = 1. + t500;
  t508 = Sin(var1[4]);
  t524 = -1.*t480*t429;
  t534 = -1.*t420*t457;
  t541 = t524 + t534;
  t494 = t429*t420;
  t496 = -1.*t480*t457;
  t497 = t494 + t496;
  t458 = -0.0695*t457;
  t571 = -1.*t429;
  t572 = 1. + t571;
  t506 = -0.0265*t502;
  t515 = -0.2375*t508;
  t517 = t506 + t515;
  t542 = -0.2375*t502;
  t544 = 0.0265*t508;
  t546 = t542 + t544;
  t587 = t480*t429;
  t593 = t420*t457;
  t595 = t587 + t593;
  t563 = t498*t497;
  t443 = 0.0265*t429;
  t467 = t443 + t458;
  t482 = -0.0695*t429;
  t486 = -0.0265*t457;
  t491 = t482 + t486;
  t601 = -1.*t595*t517;
  t602 = -1.*t497*t546;
  t605 = t498*t595;
  t606 = t497*t508;
  t607 = t605 + t606;
  t609 = 0.0265*t607;
  t611 = -1.*t595*t508;
  t614 = t563 + t611;
  t617 = -0.0225*t614;
  t573 = -0.0265*t572;
  t574 = t573 + t458;
  t580 = -0.0695*t572;
  t583 = 0.0265*t457;
  t584 = t580 + t583;
  t630 = -1.*t429*t420;
  t631 = t480*t457;
  t633 = t630 + t631;
  t578 = -1.*t480*t574;
  t586 = -1.*t420*t584;
  t618 = t578 + t586 + t601 + t602 + t609 + t617;
  t627 = t420*t574;
  t628 = -1.*t480*t584;
  t634 = -1.*t633*t517;
  t635 = -1.*t595*t546;
  t638 = -1.*t633*t508;
  t641 = t605 + t638;
  t643 = -0.0225*t641;
  t645 = t498*t633;
  t646 = t595*t508;
  t648 = t645 + t646;
  t649 = 0.0265*t648;
  t650 = t627 + t628 + t634 + t635 + t643 + t649;
  t550 = t498*t541;
  t654 = Power(t618,2);
  t655 = Power(t650,2);
  t656 = 0.00085849 + t654 + t655;
  t658 = 1/Sqrt(t656);
  t555 = -1.*t497*t508;
  t677 = 0.0265*t498;
  t678 = t677 + t515;
  t681 = -0.2375*t498;
  t682 = -0.0265*t508;
  t683 = t681 + t682;
  t704 = Cos(var1[5]);
  t706 = Sin(var1[5]);
  t709 = Cos(var1[6]);
  t711 = Sin(var1[6]);
  t705 = t480*t704;
  t707 = -1.*t420*t706;
  t708 = t705 + t707;
  t716 = -1.*t704*t420;
  t717 = -1.*t480*t706;
  t719 = t716 + t717;
  t733 = -1.*t704;
  t734 = 1. + t733;
  t721 = -0.0265*t711;
  t743 = -1.*t709;
  t744 = 1. + t743;
  t728 = t709*t719;
  t729 = -1.*t708*t711;
  t730 = t728 + t729;
  t710 = -0.0265*t709;
  t712 = -0.2375*t711;
  t713 = t710 + t712;
  t720 = 0.2375*t709;
  t722 = t720 + t721;
  t752 = t709*t708;
  t759 = t704*t420;
  t760 = t480*t706;
  t761 = t759 + t760;
  t735 = -0.0695*t734;
  t736 = -0.0265*t706;
  t737 = t735 + t736;
  t739 = -0.0265*t734;
  t740 = 0.0695*t706;
  t741 = t739 + t740;
  t745 = -0.2375*t744;
  t746 = t745 + t721;
  t748 = -0.0265*t744;
  t749 = 0.2375*t711;
  t750 = t748 + t749;
  t764 = -1.*t761*t711;
  t765 = t752 + t764;
  t738 = -1.*t480*t737;
  t742 = t420*t741;
  t747 = -1.*t708*t746;
  t751 = -1.*t719*t750;
  t753 = t719*t711;
  t754 = t752 + t753;
  t755 = -0.0225*t754;
  t756 = 0.0265*t730;
  t757 = t738 + t742 + t747 + t751 + t755 + t756;
  t771 = -1.*t420*t737;
  t772 = -1.*t480*t741;
  t773 = -1.*t761*t746;
  t774 = -1.*t708*t750;
  t775 = 0.0265*t765;
  t776 = t709*t761;
  t777 = t708*t711;
  t778 = t776 + t777;
  t779 = -0.0225*t778;
  t780 = t771 + t772 + t773 + t774 + t775 + t779;
  t783 = Power(t757,2);
  t784 = Power(t780,2);
  t785 = 0.00085849 + t783 + t784;
  t786 = 1/Sqrt(t785);
  t788 = -0.0265*t704;
  t789 = -0.0695*t706;
  t790 = t788 + t789;
  t792 = 0.0695*t704;
  t793 = t792 + t736;
  t800 = -1.*t480*t704;
  t801 = t420*t706;
  t802 = t800 + t801;
  t725 = -1.*t719*t711;
  t799 = -1.*t719*t746;
  t803 = -1.*t802*t750;
  t804 = t709*t802;
  t805 = t804 + t725;
  t806 = 0.0265*t805;
  t807 = t802*t711;
  t808 = t728 + t807;
  t809 = -0.0225*t808;
  p_output1[0]=0.5*(2.*t618*t650 + 2.*(-1.*t517*t541 + t480*t574 + t420*t584 - 1.*t546*t633 + 0.0265*(t550 + t508*t633) - 0.0225*(-1.*t508*t541 + t645))*t650)*t658*var1[9] + 0.5*(2.*(-1.*t420*t467 - 1.*t480*t491 - 1.*t497*t517 - 1.*t541*t546 - 0.0225*(t550 + t555) + 0.0265*(t508*t541 + t563))*t618 + 2.*(-1.*t467*t480 + t420*t491 + t601 + t602 + t609 + t617)*t650)*t658*var1[10] + 0.5*t658*(2.*t618*(-0.0225*(t555 - 1.*t498*t595) + 0.0265*t614 - 1.*t497*t678 - 1.*t595*t683) + 2.*t650*(-0.0225*(t611 - 1.*t498*t633) + 0.0265*t641 - 1.*t595*t678 - 1.*t633*t683))*var1[11];
  p_output1[1]=0.5*t786*(2.*t757*t780 + 2.*t757*(t420*t737 + t480*t741 + t799 + t803 + t806 + t809))*var1[9] + 0.5*t786*(2.*t780*(t747 + t751 + t755 + t756 - 1.*t420*t790 - 1.*t480*t793) + 2.*t757*(-1.*t480*t790 + t420*t793 + t799 + t803 + t806 + t809))*var1[12] + 0.5*(2.*(-1.*t708*t713 - 1.*t719*t722 + 0.0265*(-1.*t708*t709 + t725) - 0.0225*t730)*t757 + 2.*(-1.*t708*t722 - 1.*t713*t761 + 0.0265*(t729 - 1.*t709*t761) - 0.0225*t765)*t780)*t786*var1[13];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "vLegL.hh"

namespace SymFunction
{

void vLegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
