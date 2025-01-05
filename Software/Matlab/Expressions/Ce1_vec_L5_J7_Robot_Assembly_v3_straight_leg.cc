/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:42 GMT-05:00
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
  double t580;
  double t568;
  double t571;
  double t583;
  double t596;
  double t558;
  double t605;
  double t609;
  double t625;
  double t635;
  double t636;
  double t642;
  double t628;
  double t644;
  double t672;
  double t675;
  double t676;
  double t686;
  double t688;
  double t692;
  double t696;
  double t698;
  double t701;
  double t703;
  double t573;
  double t587;
  double t591;
  double t646;
  double t656;
  double t716;
  double t731;
  double t732;
  double t733;
  double t734;
  double t735;
  double t745;
  double t747;
  double t744;
  double t749;
  double t751;
  double t752;
  double t753;
  double t754;
  double t755;
  double t756;
  double t761;
  double t762;
  double t766;
  double t741;
  double t742;
  double t775;
  double t776;
  double t780;
  double t781;
  double t782;
  double t784;
  double t786;
  double t787;
  double t788;
  double t789;
  double t771;
  double t772;
  double t743;
  double t767;
  double t769;
  double t791;
  double t803;
  double t804;
  double t806;
  double t808;
  double t810;
  double t796;
  double t815;
  double t816;
  double t817;
  double t805;
  double t812;
  double t813;
  double t798;
  t580 = Cos(var1[5]);
  t568 = Cos(var1[6]);
  t571 = Sin(var1[5]);
  t583 = Sin(var1[6]);
  t596 = Cos(var1[2]);
  t558 = Sin(var1[2]);
  t605 = t580*t568;
  t609 = -1.*t571*t583;
  t625 = t605 + t609;
  t635 = -1.*t568*t571;
  t636 = -1.*t580*t583;
  t642 = t635 + t636;
  t628 = t596*t625;
  t644 = t596*t642;
  t672 = t558*t642;
  t675 = t672 + t628;
  t676 = -0.0002543413600000002*t675;
  t686 = -1.*t580*t568;
  t688 = t571*t583;
  t692 = t686 + t688;
  t696 = t558*t692;
  t698 = t644 + t696;
  t701 = -0.015373477840000005*t698;
  t703 = t676 + t701;
  t573 = t568*t571;
  t587 = t580*t583;
  t591 = t573 + t587;
  t646 = -1.*t558*t625;
  t656 = t644 + t646;
  t716 = -1.*t558*t642;
  t731 = -0.0002543413600000002*t656;
  t732 = t596*t692;
  t733 = t716 + t732;
  t734 = -0.015373477840000005*t733;
  t735 = t731 + t734;
  t745 = -1.*t568;
  t747 = 1. + t745;
  t744 = -0.0265*t571;
  t749 = -0.0265*t747;
  t751 = -0.025226*t568;
  t752 = -0.07700600000000002*t583;
  t753 = t749 + t751 + t752;
  t754 = -1.*t571*t753;
  t755 = -0.2375*t747;
  t756 = -0.314506*t568;
  t761 = -0.0012740000000000008*t583;
  t762 = t755 + t756 + t761;
  t766 = t580*t762;
  t741 = -1.*t580;
  t742 = 1. + t741;
  t775 = -0.0265*t580;
  t776 = -0.0695*t571;
  t780 = -1.*t580*t753;
  t781 = -1.*t571*t762;
  t782 = t775 + t776 + t780 + t781;
  t784 = -0.0265*t742;
  t786 = 0.0695*t571;
  t787 = t580*t753;
  t788 = t571*t762;
  t789 = t784 + t786 + t787 + t788;
  t771 = 0.0695*t580;
  t772 = t771 + t744 + t754 + t766;
  t743 = -0.0695*t742;
  t767 = t743 + t744 + t754 + t766;
  t769 = t767*t642;
  t791 = t789*t625;
  t803 = -0.07700600000000002*t568;
  t804 = t803 + t761;
  t806 = -0.0012740000000000008*t568;
  t808 = 0.07700600000000002*t583;
  t810 = t806 + t808;
  t796 = -1.*t789*t642;
  t815 = -1.*t571*t804;
  t816 = t580*t810;
  t817 = t815 + t816;
  t805 = t580*t804;
  t812 = t571*t810;
  t813 = t805 + t812;
  t798 = -1.*t767*t692;
  p_output1[0]=var2[6]*(-0.5*(-0.0002543413600000002*(-1.*t558*t591 + t628) - 0.015373477840000005*t656)*var2[2] - 0.5*t703*var2[5] - 0.5*t703*var2[6]);
  p_output1[1]=var2[6]*(-0.5*(-0.0002543413600000002*(-1.*t591*t596 + t646) - 0.015373477840000005*(-1.*t596*t625 + t716))*var2[2] - 0.5*t735*var2[5] - 0.5*t735*var2[6]);
  p_output1[2]=var2[6]*(-0.5*(-0.015373477840000005*(t769 + t591*t772 + t625*t782 + t791) - 0.0002543413600000002*(-1.*t625*t772 - 1.*t642*t782 + t796 + t798))*var2[5] - 0.5*(-0.015373477840000005*(t769 + t791 + t591*t813 + t625*t817) - 0.0002543413600000002*(t796 + t798 - 1.*t625*t813 - 1.*t642*t817))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(-0.015373477840000005*(0.0265*t568 - 0.0695*t583 + t568*t753 - 1.*t583*t762 + t583*t804 + t568*t810) - 0.0002543413600000002*(0.0695*t568 + 0.0265*t583 + t583*t753 + t568*t762 - 1.*t568*t804 + t583*t810))*Power(var2[6],2);
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

#include "Ce1_vec_L5_J7_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L5_J7_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
