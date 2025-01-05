/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:31:46 GMT-05:00
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
  double t752;
  double t736;
  double t740;
  double t755;
  double t791;
  double t744;
  double t769;
  double t776;
  double t731;
  double t836;
  double t837;
  double t787;
  double t797;
  double t803;
  double t806;
  double t808;
  double t812;
  double t846;
  double t847;
  double t849;
  double t850;
  double t840;
  double t841;
  double t842;
  double t843;
  double t816;
  double t820;
  double t821;
  double t822;
  double t825;
  double t828;
  double t876;
  double t869;
  double t870;
  double t873;
  double t874;
  double t879;
  double t884;
  double t886;
  double t903;
  double t905;
  double t906;
  double t900;
  double t901;
  double t902;
  double t907;
  double t918;
  double t919;
  double t920;
  double t921;
  double t922;
  double t912;
  double t913;
  double t914;
  double t915;
  double t916;
  double t909;
  double t910;
  double t894;
  double t896;
  double t835;
  double t844;
  double t845;
  double t851;
  double t852;
  double t857;
  double t858;
  double t859;
  double t861;
  double t862;
  double t888;
  double t889;
  double t890;
  double t887;
  double t892;
  double t893;
  double t897;
  double t898;
  double t908;
  double t911;
  double t917;
  double t923;
  double t924;
  double t926;
  double t927;
  double t928;
  double t930;
  double t931;
  double t932;
  double t934;
  double t935;
  double t936;
  double t937;
  double t938;
  double t925;
  double t929;
  double t933;
  double t939;
  double t940;
  double t957;
  double t958;
  double t959;
  double t960;
  double t961;
  double t856;
  double t864;
  double t867;
  double t948;
  double t949;
  double t950;
  double t969;
  double t970;
  double t971;
  double t814;
  double t830;
  double t831;
  double t944;
  double t945;
  double t946;
  double t965;
  double t966;
  double t967;
  t752 = Cos(var1[3]);
  t736 = Cos(var1[4]);
  t740 = Sin(var1[3]);
  t755 = Sin(var1[4]);
  t791 = Sin(var1[2]);
  t744 = -1.*t736*t740;
  t769 = -1.*t752*t755;
  t776 = t744 + t769;
  t731 = Cos(var1[2]);
  t836 = -1.*t736;
  t837 = 1. + t836;
  t787 = t731*t776;
  t797 = t752*t736;
  t803 = -1.*t740*t755;
  t806 = t797 + t803;
  t808 = t791*t806;
  t812 = t787 + t808;
  t846 = -0.0265*t837;
  t847 = -0.025229*t736;
  t849 = 0.07701400000000003*t755;
  t850 = t846 + t847 + t849;
  t840 = -0.2375*t837;
  t841 = -0.314514*t736;
  t842 = 0.0012709999999999978*t755;
  t843 = t840 + t841 + t842;
  t816 = t791*t776;
  t820 = -1.*t752*t736;
  t821 = t740*t755;
  t822 = t820 + t821;
  t825 = t731*t822;
  t828 = t816 + t825;
  t876 = t731*t806;
  t869 = t736*t740;
  t870 = t752*t755;
  t873 = t869 + t870;
  t874 = t791*t873;
  t879 = t874 + t876;
  t884 = -1.*t791*t776;
  t886 = t884 + t876;
  t903 = -0.0695*t740;
  t905 = -1.*t740*t843;
  t906 = t752*t850;
  t900 = -1.*t752;
  t901 = 1. + t900;
  t902 = -0.0265*t901;
  t907 = t902 + t903 + t905 + t906;
  t918 = -0.0695*t901;
  t919 = 0.0265*t740;
  t920 = t752*t843;
  t921 = t740*t850;
  t922 = t918 + t919 + t920 + t921;
  t912 = -0.0695*t752;
  t913 = -0.0265*t740;
  t914 = -1.*t752*t843;
  t915 = -1.*t740*t850;
  t916 = t912 + t913 + t914 + t915;
  t909 = 0.0265*t752;
  t910 = t909 + t903 + t905 + t906;
  t894 = -1.*t791*t822;
  t896 = t787 + t894;
  t835 = -0.0695*t736;
  t844 = -1.*t736*t843;
  t845 = 0.0265*t755;
  t851 = t850*t755;
  t852 = t835 + t844 + t845 + t851;
  t857 = 0.0265*t736;
  t858 = t736*t850;
  t859 = 0.0695*t755;
  t861 = t843*t755;
  t862 = t857 + t858 + t859 + t861;
  t888 = t731*t873;
  t889 = -1.*t791*t806;
  t890 = t888 + t889;
  t887 = 0.19964*t886*t879;
  t892 = 0.19964*t890*t812;
  t893 = 0.19964*t886*t828;
  t897 = 0.19964*t812*t896;
  t898 = t887 + t892 + t893 + t897;
  t908 = -1.*t907*t776;
  t911 = -1.*t910*t873;
  t917 = -1.*t916*t806;
  t923 = -1.*t922*t806;
  t924 = t908 + t911 + t917 + t923;
  t926 = t907*t776;
  t927 = t922*t806;
  t928 = t926 + t927;
  t930 = -1.*t922*t873;
  t931 = -1.*t907*t806;
  t932 = t930 + t931;
  t934 = t916*t776;
  t935 = t922*t776;
  t936 = t910*t806;
  t937 = t907*t822;
  t938 = t934 + t935 + t936 + t937;
  t925 = 0.19964*t812*t924;
  t929 = 0.19964*t812*t928;
  t933 = 0.19964*t932*t828;
  t939 = 0.19964*t879*t938;
  t940 = t925 + t929 + t933 + t939;
  t957 = 0.19964*t886*t924;
  t958 = 0.19964*t886*t928;
  t959 = 0.19964*t932*t896;
  t960 = 0.19964*t890*t938;
  t961 = t957 + t958 + t959 + t960;
  t856 = 0.19964*t852*t812;
  t864 = 0.19964*t862*t828;
  t867 = t856 + t864;
  t948 = 0.19964*t852*t886;
  t949 = 0.19964*t862*t896;
  t950 = t948 + t949;
  t969 = 0.19964*t862*t924;
  t970 = 0.19964*t852*t938;
  t971 = t969 + t970;
  t814 = 0.015375074960000006*t812;
  t830 = 0.0002537424399999996*t828;
  t831 = t814 + t830;
  t944 = 0.015375074960000006*t886;
  t945 = 0.0002537424399999996*t896;
  t946 = t944 + t945;
  t965 = 0.0002537424399999996*t924;
  t966 = 0.015375074960000006*t938;
  t967 = t965 + t966;
  p_output1[0]=var2[3]*(-0.5*(0.39928*t812*t828 + 0.39928*t812*t879)*var2[0] - 0.5*t898*var2[1] - 0.5*t940*var2[2] - 0.5*t867*var2[3] - 0.5*t831*var2[4]);
  p_output1[1]=var2[3]*(-0.5*t898*var2[0] - 0.5*(0.39928*t886*t890 + 0.39928*t886*t896)*var2[1] - 0.5*t961*var2[2] - 0.5*t950*var2[3] - 0.5*t946*var2[4]);
  p_output1[2]=var2[3]*(-0.5*t940*var2[0] - 0.5*t961*var2[1] - 0.5*(0.39928*t924*t932 + 0.39928*t928*t938)*var2[2] - 0.5*t971*var2[3] - 0.5*t967*var2[4]);
  p_output1[3]=(-0.5*t867*var2[0] - 0.5*t950*var2[1] - 0.5*t971*var2[2])*var2[3];
  p_output1[4]=(-0.5*t831*var2[0] - 0.5*t946*var2[1] - 0.5*t967*var2[2])*var2[3];
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

#include "Ce2_vec_L3_J4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L3_J4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
