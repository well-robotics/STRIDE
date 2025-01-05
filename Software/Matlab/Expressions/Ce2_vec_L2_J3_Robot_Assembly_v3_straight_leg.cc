/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:31:24 GMT-05:00
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
  double t632;
  double t587;
  double t594;
  double t644;
  double t672;
  double t675;
  double t676;
  double t605;
  double t656;
  double t669;
  double t701;
  double t705;
  double t706;
  double t738;
  double t739;
  double t752;
  double t755;
  double t756;
  double t761;
  double t740;
  double t743;
  double t744;
  double t749;
  double t728;
  double t729;
  double t730;
  double t714;
  double t724;
  double t725;
  double t727;
  double t731;
  double t732;
  double t733;
  double t734;
  double t751;
  double t766;
  double t769;
  double t774;
  double t775;
  double t776;
  double t771;
  double t783;
  double t784;
  double t802;
  double t803;
  double t805;
  double t791;
  double t793;
  double t794;
  double t671;
  double t688;
  double t696;
  t632 = Cos(var1[2]);
  t587 = Cos(var1[3]);
  t594 = Sin(var1[2]);
  t644 = Sin(var1[3]);
  t672 = t632*t587;
  t675 = t594*t644;
  t676 = t672 + t675;
  t605 = -1.*t587*t594;
  t656 = t632*t644;
  t669 = t605 + t656;
  t701 = t587*t594;
  t705 = -1.*t632*t644;
  t706 = t701 + t705;
  t738 = -1.*t587;
  t739 = 1. + t738;
  t752 = -0.0695*t739;
  t755 = -0.15232*t587;
  t756 = 0.0011329999999999986*t644;
  t761 = t752 + t755 + t756;
  t740 = -0.0265*t739;
  t743 = -0.025367*t587;
  t744 = 0.08282*t644;
  t749 = t740 + t743 + t744;
  t728 = -1.*t632*t587;
  t729 = -1.*t594*t644;
  t730 = t728 + t729;
  t714 = 1.38102*t669*t676;
  t724 = 0.69051*t706*t669;
  t725 = Power(t669,2);
  t727 = 0.69051*t725;
  t731 = 0.69051*t730*t676;
  t732 = Power(t676,2);
  t733 = 0.69051*t732;
  t734 = t724 + t727 + t731 + t733;
  t751 = -1.*t587*t749;
  t766 = -1.*t761*t644;
  t769 = t751 + t766;
  t774 = t587*t761;
  t775 = -1.*t749*t644;
  t776 = t774 + t775;
  t771 = 0.69051*t676*t769;
  t783 = 0.69051*t669*t776;
  t784 = t771 + t783;
  t802 = 0.69051*t669*t769;
  t803 = 0.69051*t730*t776;
  t805 = t802 + t803;
  t791 = 0.0007823478299999989*t669;
  t793 = 0.0571880382*t730;
  t794 = t791 + t793;
  t671 = 0.0571880382*t669;
  t688 = 0.0007823478299999989*t676;
  t696 = t671 + t688;
  p_output1[0]=var2[2]*(-0.5*(1.38102*t676*t706 + t714)*var2[0] - 0.5*t734*var2[1] - 0.5*t784*var2[2] - 0.5*t696*var2[3]);
  p_output1[1]=var2[2]*(-0.5*t734*var2[0] - 0.5*(t714 + 1.38102*t669*t730)*var2[1] - 0.5*t805*var2[2] - 0.5*t794*var2[3]);
  p_output1[2]=(-0.5*t784*var2[0] - 0.5*t805*var2[1])*var2[2];
  p_output1[3]=(-0.5*t696*var2[0] - 0.5*t794*var2[1])*var2[2];
  p_output1[4]=0;
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

#include "Ce2_vec_L2_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L2_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
