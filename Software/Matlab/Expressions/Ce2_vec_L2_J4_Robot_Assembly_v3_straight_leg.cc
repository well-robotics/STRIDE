/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:31:27 GMT-05:00
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
  double t672;
  double t656;
  double t669;
  double t675;
  double t671;
  double t676;
  double t688;
  double t701;
  double t706;
  double t707;
  double t725;
  double t727;
  double t730;
  double t783;
  double t786;
  double t787;
  double t788;
  double t791;
  double t793;
  double t795;
  double t797;
  double t799;
  double t800;
  double t801;
  double t802;
  double t814;
  double t816;
  double t818;
  double t808;
  double t810;
  double t731;
  double t743;
  double t744;
  double t752;
  double t736;
  double t740;
  double t755;
  double t756;
  double t769;
  double t771;
  double t774;
  double t803;
  double t812;
  double t819;
  double t820;
  double t822;
  double t823;
  double t824;
  double t826;
  double t827;
  double t828;
  double t829;
  double t830;
  double t806;
  double t821;
  double t825;
  double t831;
  double t832;
  double t844;
  double t845;
  double t846;
  double t847;
  double t848;
  double t698;
  double t714;
  double t718;
  double t836;
  double t837;
  double t838;
  double t852;
  double t853;
  double t854;
  t672 = Cos(var1[2]);
  t656 = Cos(var1[3]);
  t669 = Sin(var1[2]);
  t675 = Sin(var1[3]);
  t671 = t656*t669;
  t676 = -1.*t672*t675;
  t688 = t671 + t676;
  t701 = -1.*t672*t656;
  t706 = -1.*t669*t675;
  t707 = t701 + t706;
  t725 = t672*t656;
  t727 = t669*t675;
  t730 = t725 + t727;
  t783 = -1.*t656;
  t786 = 1. + t783;
  t787 = -0.0265*t786;
  t788 = -0.025367*t656;
  t791 = 0.08282*t675;
  t793 = t787 + t788 + t791;
  t795 = -1.*t656*t793;
  t797 = -0.0695*t786;
  t799 = -0.15232*t656;
  t800 = 0.0011329999999999986*t675;
  t801 = t797 + t799 + t800;
  t802 = -1.*t801*t675;
  t814 = 0.08282*t656;
  t816 = -0.0011329999999999986*t675;
  t818 = t814 + t816;
  t808 = 0.0011329999999999986*t656;
  t810 = t808 + t791;
  t731 = 1.38102*t688*t730;
  t743 = -1.*t656*t669;
  t744 = t672*t675;
  t752 = t743 + t744;
  t736 = Power(t688,2);
  t740 = 0.69051*t736;
  t755 = 0.69051*t688*t752;
  t756 = 0.69051*t707*t730;
  t769 = Power(t730,2);
  t771 = 0.69051*t769;
  t774 = t740 + t755 + t756 + t771;
  t803 = t795 + t802;
  t812 = t656*t810;
  t819 = -1.*t818*t675;
  t820 = t795 + t812 + t819 + t802;
  t822 = t656*t801;
  t823 = -1.*t793*t675;
  t824 = t822 + t823;
  t826 = -1.*t656*t818;
  t827 = -1.*t656*t801;
  t828 = t793*t675;
  t829 = -1.*t810*t675;
  t830 = t826 + t827 + t828 + t829;
  t806 = 0.69051*t707*t803;
  t821 = 0.69051*t730*t820;
  t825 = 0.69051*t688*t824;
  t831 = 0.69051*t688*t830;
  t832 = t806 + t821 + t825 + t831;
  t844 = 0.69051*t688*t803;
  t845 = 0.69051*t752*t820;
  t846 = 0.69051*t730*t824;
  t847 = 0.69051*t730*t830;
  t848 = t844 + t845 + t846 + t847;
  t698 = 0.0571880382*t688;
  t714 = 0.0007823478299999989*t707;
  t718 = t698 + t714;
  t836 = 0.0007823478299999989*t688;
  t837 = 0.0571880382*t730;
  t838 = t836 + t837;
  t852 = 0.0571880382*t820;
  t853 = 0.0007823478299999989*t830;
  t854 = t852 + t853;
  p_output1[0]=var2[3]*(-0.5*(1.38102*t688*t707 + t731)*var2[0] - 0.5*t774*var2[1] - 0.5*t832*var2[2] - 0.5*t718*var2[3]);
  p_output1[1]=var2[3]*(-0.5*t774*var2[0] - 0.5*(t731 + 1.38102*t730*t752)*var2[1] - 0.5*t848*var2[2] - 0.5*t838*var2[3]);
  p_output1[2]=var2[3]*(-0.5*t832*var2[0] - 0.5*t848*var2[1] - 0.5*(1.38102*t820*t824 + 1.38102*t803*t830)*var2[2] - 0.5*t854*var2[3]);
  p_output1[3]=(-0.5*t718*var2[0] - 0.5*t838*var2[1] - 0.5*t854*var2[2])*var2[3];
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

#include "Ce2_vec_L2_J4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L2_J4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
