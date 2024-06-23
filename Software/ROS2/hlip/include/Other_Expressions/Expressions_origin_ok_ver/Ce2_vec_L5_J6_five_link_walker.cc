/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:47 GMT-05:00
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
  double t793;
  double t765;
  double t775;
  double t794;
  double t810;
  double t783;
  double t803;
  double t806;
  double t755;
  double t807;
  double t816;
  double t820;
  double t829;
  double t830;
  double t831;
  double t858;
  double t860;
  double t861;
  double t863;
  double t869;
  double t873;
  double t852;
  double t854;
  double t855;
  double t834;
  double t839;
  double t840;
  double t842;
  double t843;
  double t844;
  double t891;
  double t893;
  double t899;
  double t900;
  double t902;
  double t905;
  double t917;
  double t918;
  double t931;
  double t932;
  double t933;
  double t935;
  double t937;
  double t938;
  double t944;
  double t945;
  double t946;
  double t921;
  double t922;
  double t923;
  double t857;
  double t874;
  double t877;
  double t879;
  double t881;
  double t883;
  double t912;
  double t914;
  double t916;
  double t920;
  double t925;
  double t926;
  double t927;
  double t934;
  double t939;
  double t940;
  double t942;
  double t943;
  double t947;
  double t949;
  double t950;
  double t952;
  double t953;
  double t954;
  double t956;
  double t957;
  double t958;
  double t959;
  double t960;
  double t979;
  double t980;
  double t981;
  double t982;
  double t983;
  double t941;
  double t951;
  double t955;
  double t961;
  double t962;
  double t970;
  double t971;
  double t972;
  double t878;
  double t885;
  double t888;
  double t991;
  double t992;
  double t993;
  double t966;
  double t967;
  double t968;
  double t832;
  double t846;
  double t849;
  double t987;
  double t988;
  double t989;
  t793 = Cos(var1[5]);
  t765 = Cos(var1[6]);
  t775 = Sin(var1[5]);
  t794 = Sin(var1[6]);
  t810 = Cos(var1[2]);
  t783 = -1.*t765*t775;
  t803 = -1.*t793*t794;
  t806 = t783 + t803;
  t755 = Sin(var1[2]);
  t807 = t755*t806;
  t816 = t793*t765;
  t820 = -1.*t775*t794;
  t829 = t816 + t820;
  t830 = t810*t829;
  t831 = t807 + t830;
  t858 = -1.*t765;
  t860 = 1. + t858;
  t861 = -0.16*t860;
  t863 = -0.167368*t765;
  t869 = 0.022659*t794;
  t873 = t861 + t863 + t869;
  t852 = -0.022659*t765;
  t854 = -0.007367999999999986*t794;
  t855 = t852 + t854;
  t834 = t810*t806;
  t839 = -1.*t793*t765;
  t840 = t775*t794;
  t842 = t839 + t840;
  t843 = t755*t842;
  t844 = t834 + t843;
  t891 = t765*t775;
  t893 = t793*t794;
  t899 = t891 + t893;
  t900 = t810*t899;
  t902 = t755*t829;
  t905 = t900 + t902;
  t917 = -1.*t755*t829;
  t918 = t834 + t917;
  t931 = -1.*t775*t855;
  t932 = t793*t873;
  t933 = t931 + t932;
  t935 = t793*t855;
  t937 = t775*t873;
  t938 = t935 + t937;
  t944 = -1.*t793*t855;
  t945 = -1.*t775*t873;
  t946 = t944 + t945;
  t921 = -1.*t755*t806;
  t922 = t810*t842;
  t923 = t921 + t922;
  t857 = -1.*t765*t855;
  t874 = t873*t794;
  t877 = t857 + t874;
  t879 = t765*t873;
  t881 = t855*t794;
  t883 = t879 + t881;
  t912 = -1.*t755*t899;
  t914 = t912 + t830;
  t916 = 0.14994*t831*t914;
  t920 = 0.14994*t918*t905;
  t925 = 0.14994*t831*t923;
  t926 = 0.14994*t918*t844;
  t927 = t916 + t920 + t925 + t926;
  t934 = -1.*t933*t806;
  t939 = -1.*t938*t829;
  t940 = t934 + t939;
  t942 = t933*t806;
  t943 = t933*t899;
  t947 = t946*t829;
  t949 = t938*t829;
  t950 = t942 + t943 + t947 + t949;
  t952 = t938*t899;
  t953 = t933*t829;
  t954 = t952 + t953;
  t956 = -1.*t946*t806;
  t957 = -1.*t938*t806;
  t958 = -1.*t933*t829;
  t959 = -1.*t933*t842;
  t960 = t956 + t957 + t958 + t959;
  t979 = 0.14994*t918*t940;
  t980 = 0.14994*t918*t950;
  t981 = 0.14994*t954*t923;
  t982 = 0.14994*t914*t960;
  t983 = t979 + t980 + t981 + t982;
  t941 = 0.14994*t831*t940;
  t951 = 0.14994*t831*t950;
  t955 = 0.14994*t954*t844;
  t961 = 0.14994*t905*t960;
  t962 = t941 + t951 + t955 + t961;
  t970 = 0.14994*t877*t918;
  t971 = 0.14994*t883*t923;
  t972 = t970 + t971;
  t878 = 0.14994*t877*t831;
  t885 = 0.14994*t883*t844;
  t888 = t878 + t885;
  t991 = 0.14994*t883*t950;
  t992 = 0.14994*t877*t960;
  t993 = t991 + t992;
  t966 = 0.0033974904599999994*t918;
  t967 = -0.0011047579199999977*t923;
  t968 = t966 + t967;
  t832 = 0.0033974904599999994*t831;
  t846 = -0.0011047579199999977*t844;
  t849 = t832 + t846;
  t987 = -0.0011047579199999977*t950;
  t988 = 0.0033974904599999994*t960;
  t989 = t987 + t988;
  p_output1[0]=var2[5]*(-0.5*(0.29988*t831*t844 + 0.29988*t831*t905)*var2[0] - 0.5*t927*var2[1] - 0.5*t962*var2[2] - 0.5*t888*var2[5] - 0.5*t849*var2[6]);
  p_output1[1]=var2[5]*(-0.5*t927*var2[0] - 0.5*(0.29988*t914*t918 + 0.29988*t918*t923)*var2[1] - 0.5*t983*var2[2] - 0.5*t972*var2[5] - 0.5*t968*var2[6]);
  p_output1[2]=var2[5]*(-0.5*t962*var2[0] - 0.5*t983*var2[1] - 0.5*(0.29988*t950*t954 + 0.29988*t940*t960)*var2[2] - 0.5*t993*var2[5] - 0.5*t989*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t888*var2[0] - 0.5*t972*var2[1] - 0.5*t993*var2[2])*var2[5];
  p_output1[6]=(-0.5*t849*var2[0] - 0.5*t968*var2[1] - 0.5*t989*var2[2])*var2[5];
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

#include "Ce2_vec_L5_J6_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L5_J6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
