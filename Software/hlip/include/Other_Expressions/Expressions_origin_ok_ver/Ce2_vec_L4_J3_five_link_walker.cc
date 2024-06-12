/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:26 GMT-05:00
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
  double t666;
  double t651;
  double t652;
  double t680;
  double t689;
  double t648;
  double t693;
  double t695;
  double t696;
  double t662;
  double t682;
  double t684;
  double t688;
  double t699;
  double t703;
  double t733;
  double t734;
  double t735;
  double t736;
  double t738;
  double t739;
  double t728;
  double t729;
  double t730;
  double t706;
  double t709;
  double t713;
  double t716;
  double t717;
  double t720;
  double t753;
  double t754;
  double t756;
  double t757;
  double t759;
  double t782;
  double t783;
  double t784;
  double t778;
  double t779;
  double t780;
  double t763;
  double t764;
  double t765;
  double t772;
  double t773;
  double t746;
  double t747;
  double t748;
  double t732;
  double t741;
  double t742;
  double t766;
  double t767;
  double t769;
  double t770;
  double t771;
  double t775;
  double t776;
  double t781;
  double t785;
  double t786;
  double t788;
  double t789;
  double t790;
  double t787;
  double t791;
  double t792;
  double t809;
  double t810;
  double t811;
  double t744;
  double t749;
  double t750;
  double t800;
  double t801;
  double t802;
  double t704;
  double t722;
  double t723;
  double t796;
  double t797;
  double t798;
  t666 = Cos(var1[3]);
  t651 = Cos(var1[4]);
  t652 = Sin(var1[3]);
  t680 = Sin(var1[4]);
  t689 = Cos(var1[2]);
  t648 = Sin(var1[2]);
  t693 = t666*t651;
  t695 = -1.*t652*t680;
  t696 = t693 + t695;
  t662 = t651*t652;
  t682 = t666*t680;
  t684 = t662 + t682;
  t688 = -1.*t648*t684;
  t699 = t689*t696;
  t703 = t688 + t699;
  t733 = -1.*t651;
  t734 = 1. + t733;
  t735 = -0.16*t734;
  t736 = -0.167371*t651;
  t738 = 0.022663*t680;
  t739 = t735 + t736 + t738;
  t728 = -0.022663*t651;
  t729 = -0.007370999999999989*t680;
  t730 = t728 + t729;
  t706 = -1.*t651*t652;
  t709 = -1.*t666*t680;
  t713 = t706 + t709;
  t716 = t689*t713;
  t717 = -1.*t648*t696;
  t720 = t716 + t717;
  t753 = t648*t713;
  t754 = t753 + t699;
  t756 = t689*t684;
  t757 = t648*t696;
  t759 = t756 + t757;
  t782 = -1.*t652*t730;
  t783 = t666*t739;
  t784 = t782 + t783;
  t778 = t666*t730;
  t779 = t652*t739;
  t780 = t778 + t779;
  t763 = -1.*t648*t713;
  t764 = -1.*t689*t696;
  t765 = t763 + t764;
  t772 = -1.*t689*t684;
  t773 = t772 + t717;
  t746 = t651*t739;
  t747 = t730*t680;
  t748 = t746 + t747;
  t732 = -1.*t651*t730;
  t741 = t739*t680;
  t742 = t732 + t741;
  t766 = 0.14994*t765*t754;
  t767 = Power(t703,2);
  t769 = 0.14994*t767;
  t770 = Power(t720,2);
  t771 = 0.14994*t770;
  t775 = 0.14994*t773*t759;
  t776 = t766 + t769 + t771 + t775;
  t781 = t780*t684;
  t785 = t784*t696;
  t786 = t781 + t785;
  t788 = -1.*t784*t713;
  t789 = -1.*t780*t696;
  t790 = t788 + t789;
  t787 = 0.14994*t720*t786;
  t791 = 0.14994*t703*t790;
  t792 = t787 + t791;
  t809 = 0.14994*t765*t786;
  t810 = 0.14994*t773*t790;
  t811 = t809 + t810;
  t744 = 0.14994*t742*t703;
  t749 = 0.14994*t748*t720;
  t750 = t744 + t749;
  t800 = 0.14994*t748*t765;
  t801 = 0.14994*t742*t773;
  t802 = t800 + t801;
  t704 = 0.0033980902199999994*t703;
  t722 = -0.0011052077399999983*t720;
  t723 = t704 + t722;
  t796 = -0.0011052077399999983*t765;
  t797 = 0.0033980902199999994*t773;
  t798 = t796 + t797;
  p_output1[0]=var2[2]*(-0.5*(0.29988*t720*t754 + 0.29988*t703*t759)*var2[0] - 0.5*t776*var2[1] - 0.5*t792*var2[2] - 0.5*t750*var2[3] - 0.5*t723*var2[4]);
  p_output1[1]=var2[2]*(-0.5*t776*var2[0] - 0.5*(0.29988*t720*t765 + 0.29988*t703*t773)*var2[1] - 0.5*t811*var2[2] - 0.5*t802*var2[3] - 0.5*t798*var2[4]);
  p_output1[2]=(-0.5*t792*var2[0] - 0.5*t811*var2[1])*var2[2];
  p_output1[3]=(-0.5*t750*var2[0] - 0.5*t802*var2[1])*var2[2];
  p_output1[4]=(-0.5*t723*var2[0] - 0.5*t798*var2[1])*var2[2];
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

#include "Ce2_vec_L4_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L4_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
