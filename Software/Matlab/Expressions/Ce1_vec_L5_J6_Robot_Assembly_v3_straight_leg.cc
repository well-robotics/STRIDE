/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:40 GMT-05:00
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
  double t555;
  double t572;
  double t563;
  double t567;
  double t606;
  double t600;
  double t583;
  double t587;
  double t591;
  double t592;
  double t568;
  double t571;
  double t573;
  double t578;
  double t623;
  double t598;
  double t625;
  double t626;
  double t627;
  double t558;
  double t579;
  double t580;
  double t594;
  double t596;
  double t653;
  double t654;
  double t655;
  double t628;
  double t635;
  double t636;
  double t642;
  double t644;
  double t646;
  double t656;
  double t675;
  double t676;
  double t683;
  double t718;
  double t719;
  double t705;
  double t706;
  double t707;
  double t605;
  double t609;
  double t613;
  double t684;
  double t686;
  double t687;
  double t688;
  double t692;
  double t694;
  double t659;
  double t669;
  double t737;
  double t726;
  double t727;
  double t728;
  double t614;
  double t632;
  double t746;
  double t701;
  double t703;
  double t714;
  double t716;
  double t717;
  double t720;
  double t724;
  double t747;
  double t749;
  double t750;
  double t767;
  double t768;
  double t769;
  double t763;
  double t765;
  double t775;
  double t776;
  double t777;
  double t778;
  double t779;
  double t781;
  double t782;
  double t783;
  double t784;
  double t785;
  double t772;
  double t773;
  double t766;
  double t770;
  double t771;
  double t786;
  double t790;
  double t809;
  double t810;
  double t811;
  double t805;
  double t806;
  double t807;
  double t792;
  t555 = Cos(var1[6]);
  t572 = Sin(var1[6]);
  t563 = -1.*t555;
  t567 = 1. + t563;
  t606 = Cos(var1[5]);
  t600 = Sin(var1[5]);
  t583 = -0.2375*t567;
  t587 = -0.314506*t555;
  t591 = -0.0012740000000000008*t572;
  t592 = t583 + t587 + t591;
  t568 = -0.0265*t567;
  t571 = -0.025226*t555;
  t573 = -0.07700600000000002*t572;
  t578 = t568 + t571 + t573;
  t623 = Cos(var1[2]);
  t598 = Sin(var1[2]);
  t625 = t606*t555;
  t626 = -1.*t600*t572;
  t627 = t625 + t626;
  t558 = -0.0265*t555;
  t579 = -1.*t555*t578;
  t580 = 0.0695*t572;
  t594 = t592*t572;
  t596 = t558 + t579 + t580 + t594;
  t653 = -1.*t555*t600;
  t654 = -1.*t606*t572;
  t655 = t653 + t654;
  t628 = t623*t627;
  t635 = 0.0695*t555;
  t636 = t555*t592;
  t642 = 0.0265*t572;
  t644 = t578*t572;
  t646 = t635 + t636 + t642 + t644;
  t656 = t623*t655;
  t675 = t598*t655;
  t676 = t675 + t628;
  t683 = 0.19964*t596*t676;
  t718 = -0.07700600000000002*t555;
  t719 = t718 + t591;
  t705 = -0.0012740000000000008*t555;
  t706 = 0.07700600000000002*t572;
  t707 = t705 + t706;
  t605 = t555*t600;
  t609 = t606*t572;
  t613 = t605 + t609;
  t684 = -1.*t606*t555;
  t686 = t600*t572;
  t687 = t684 + t686;
  t688 = t598*t687;
  t692 = t656 + t688;
  t694 = 0.19964*t646*t692;
  t659 = -1.*t598*t627;
  t669 = t656 + t659;
  t737 = -1.*t598*t655;
  t726 = -1.*t555*t719;
  t727 = t707*t572;
  t728 = t635 + t636 + t726 + t642 + t644 + t727;
  t614 = -1.*t598*t613;
  t632 = t614 + t628;
  t746 = 0.19964*t596*t669;
  t701 = 0.0265*t555;
  t703 = t555*t578;
  t714 = t555*t707;
  t716 = -0.0695*t572;
  t717 = -1.*t592*t572;
  t720 = t719*t572;
  t724 = t701 + t703 + t714 + t716 + t717 + t720;
  t747 = t623*t687;
  t749 = t737 + t747;
  t750 = 0.19964*t646*t749;
  t767 = -0.0265*t600;
  t768 = -1.*t600*t578;
  t769 = t606*t592;
  t763 = -1.*t606;
  t765 = 1. + t763;
  t775 = -0.0265*t606;
  t776 = -0.0695*t600;
  t777 = -1.*t606*t578;
  t778 = -1.*t600*t592;
  t779 = t775 + t776 + t777 + t778;
  t781 = -0.0265*t765;
  t782 = 0.0695*t600;
  t783 = t606*t578;
  t784 = t600*t592;
  t785 = t781 + t782 + t783 + t784;
  t772 = 0.0695*t606;
  t773 = t772 + t767 + t768 + t769;
  t766 = -0.0695*t765;
  t770 = t766 + t767 + t768 + t769;
  t771 = t770*t655;
  t786 = t785*t627;
  t790 = -1.*t785*t655;
  t809 = -1.*t600*t719;
  t810 = t606*t707;
  t811 = t809 + t810;
  t805 = t606*t719;
  t806 = t600*t707;
  t807 = t805 + t806;
  t792 = -1.*t770*t687;
  p_output1[0]=var2[5]*(-0.5*(0.19964*t596*t632 + 0.19964*t646*t669)*var2[2] - 0.5*(t683 + t694)*var2[5] - 0.5*(t683 + t694 + 0.19964*t676*t724 + 0.19964*(t613*t623 + t598*t627)*t728)*var2[6]);
  p_output1[1]=var2[5]*(-0.5*(0.19964*t596*(-1.*t613*t623 + t659) + 0.19964*t646*(-1.*t623*t627 + t737))*var2[2] - 0.5*(t746 + t750)*var2[5] - 0.5*(0.19964*t669*t724 + 0.19964*t632*t728 + t746 + t750)*var2[6]);
  p_output1[2]=var2[5]*(-0.5*(0.19964*t646*(t771 + t613*t773 + t627*t779 + t786) + 0.19964*t596*(-1.*t627*t773 - 1.*t655*t779 + t790 + t792))*var2[5] - 0.5*(0.19964*t724*(t627*t770 + t613*t785) + 0.19964*t728*(-1.*t655*t770 - 1.*t627*t785) + 0.19964*t646*(t771 + t786 + t613*t807 + t627*t811) + 0.19964*t596*(t790 + t792 - 1.*t627*t807 - 1.*t655*t811))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(0.39928*t646*t724 + 0.39928*t596*t728)*var2[5]*var2[6];
  p_output1[6]=-0.5*(-0.015373477840000005*t724 - 0.0002543413600000002*t728)*var2[5]*var2[6];
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

#include "Ce1_vec_L5_J6_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L5_J6_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
