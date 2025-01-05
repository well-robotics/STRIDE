/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:31:43 GMT-05:00
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
  double t714;
  double t688;
  double t698;
  double t720;
  double t730;
  double t671;
  double t731;
  double t732;
  double t733;
  double t799;
  double t800;
  double t707;
  double t724;
  double t725;
  double t727;
  double t736;
  double t740;
  double t820;
  double t821;
  double t822;
  double t824;
  double t803;
  double t806;
  double t808;
  double t812;
  double t752;
  double t755;
  double t756;
  double t769;
  double t771;
  double t776;
  double t841;
  double t842;
  double t844;
  double t845;
  double t846;
  double t865;
  double t866;
  double t873;
  double t874;
  double t875;
  double t876;
  double t877;
  double t867;
  double t868;
  double t869;
  double t870;
  double t871;
  double t853;
  double t855;
  double t856;
  double t860;
  double t861;
  double t831;
  double t833;
  double t834;
  double t835;
  double t836;
  double t797;
  double t814;
  double t816;
  double t825;
  double t828;
  double t851;
  double t852;
  double t857;
  double t858;
  double t859;
  double t862;
  double t863;
  double t872;
  double t878;
  double t879;
  double t881;
  double t882;
  double t883;
  double t902;
  double t903;
  double t904;
  double t880;
  double t884;
  double t885;
  double t893;
  double t894;
  double t895;
  double t830;
  double t837;
  double t839;
  double t889;
  double t890;
  double t891;
  double t744;
  double t787;
  double t788;
  t714 = Cos(var1[3]);
  t688 = Cos(var1[4]);
  t698 = Sin(var1[3]);
  t720 = Sin(var1[4]);
  t730 = Cos(var1[2]);
  t671 = Sin(var1[2]);
  t731 = t714*t688;
  t732 = -1.*t698*t720;
  t733 = t731 + t732;
  t799 = -1.*t688;
  t800 = 1. + t799;
  t707 = -1.*t688*t698;
  t724 = -1.*t714*t720;
  t725 = t707 + t724;
  t727 = -1.*t671*t725;
  t736 = t730*t733;
  t740 = t727 + t736;
  t820 = -0.2375*t800;
  t821 = -0.314514*t688;
  t822 = 0.0012709999999999978*t720;
  t824 = t820 + t821 + t822;
  t803 = -0.0265*t800;
  t806 = -0.025229*t688;
  t808 = 0.07701400000000003*t720;
  t812 = t803 + t806 + t808;
  t752 = t688*t698;
  t755 = t714*t720;
  t756 = t752 + t755;
  t769 = t730*t756;
  t771 = -1.*t671*t733;
  t776 = t769 + t771;
  t841 = t671*t756;
  t842 = t841 + t736;
  t844 = t730*t725;
  t845 = t671*t733;
  t846 = t844 + t845;
  t865 = -1.*t714;
  t866 = 1. + t865;
  t873 = -0.0265*t866;
  t874 = -0.0695*t698;
  t875 = -1.*t698*t824;
  t876 = t714*t812;
  t877 = t873 + t874 + t875 + t876;
  t867 = -0.0695*t866;
  t868 = 0.0265*t698;
  t869 = t714*t824;
  t870 = t698*t812;
  t871 = t867 + t868 + t869 + t870;
  t853 = -1.*t671*t756;
  t855 = -1.*t730*t733;
  t856 = t853 + t855;
  t860 = -1.*t730*t725;
  t861 = t860 + t771;
  t831 = -0.0695*t688;
  t833 = -1.*t688*t824;
  t834 = 0.0265*t720;
  t835 = t812*t720;
  t836 = t831 + t833 + t834 + t835;
  t797 = 0.0265*t688;
  t814 = t688*t812;
  t816 = 0.0695*t720;
  t825 = t824*t720;
  t828 = t797 + t814 + t816 + t825;
  t851 = Power(t740,2);
  t852 = 0.19964*t851;
  t857 = 0.19964*t856*t842;
  t858 = Power(t776,2);
  t859 = 0.19964*t858;
  t862 = 0.19964*t861*t846;
  t863 = t852 + t857 + t859 + t862;
  t872 = -1.*t871*t756;
  t878 = -1.*t877*t733;
  t879 = t872 + t878;
  t881 = t877*t725;
  t882 = t871*t733;
  t883 = t881 + t882;
  t902 = 0.19964*t861*t879;
  t903 = 0.19964*t856*t883;
  t904 = t902 + t903;
  t880 = 0.19964*t740*t879;
  t884 = 0.19964*t776*t883;
  t885 = t880 + t884;
  t893 = 0.19964*t836*t856;
  t894 = 0.19964*t828*t861;
  t895 = t893 + t894;
  t830 = 0.19964*t828*t740;
  t837 = 0.19964*t836*t776;
  t839 = t830 + t837;
  t889 = 0.015375074960000006*t856;
  t890 = 0.0002537424399999996*t861;
  t891 = t889 + t890;
  t744 = 0.0002537424399999996*t740;
  t787 = 0.015375074960000006*t776;
  t788 = t744 + t787;
  p_output1[0]=var2[2]*(-0.5*(0.39928*t776*t842 + 0.39928*t740*t846)*var2[0] - 0.5*t863*var2[1] - 0.5*t885*var2[2] - 0.5*t839*var2[3] - 0.5*t788*var2[4]);
  p_output1[1]=var2[2]*(-0.5*t863*var2[0] - 0.5*(0.39928*t776*t856 + 0.39928*t740*t861)*var2[1] - 0.5*t904*var2[2] - 0.5*t895*var2[3] - 0.5*t891*var2[4]);
  p_output1[2]=(-0.5*t885*var2[0] - 0.5*t904*var2[1])*var2[2];
  p_output1[3]=(-0.5*t839*var2[0] - 0.5*t895*var2[1])*var2[2];
  p_output1[4]=(-0.5*t788*var2[0] - 0.5*t891*var2[1])*var2[2];
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

#include "Ce2_vec_L3_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L3_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
