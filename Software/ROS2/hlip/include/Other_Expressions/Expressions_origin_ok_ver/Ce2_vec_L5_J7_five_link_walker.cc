/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:49 GMT-05:00
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
  double t831;
  double t810;
  double t816;
  double t832;
  double t844;
  double t830;
  double t834;
  double t840;
  double t807;
  double t843;
  double t846;
  double t850;
  double t852;
  double t854;
  double t861;
  double t891;
  double t893;
  double t900;
  double t905;
  double t906;
  double t909;
  double t910;
  double t911;
  double t914;
  double t935;
  double t937;
  double t923;
  double t925;
  double t926;
  double t943;
  double t949;
  double t869;
  double t874;
  double t877;
  double t878;
  double t879;
  double t881;
  double t952;
  double t953;
  double t954;
  double t955;
  double t960;
  double t961;
  double t977;
  double t978;
  double t988;
  double t990;
  double t991;
  double t994;
  double t995;
  double t996;
  double t1006;
  double t1007;
  double t1008;
  double t1001;
  double t1002;
  double t1003;
  double t980;
  double t981;
  double t982;
  double t947;
  double t950;
  double t951;
  double t974;
  double t975;
  double t902;
  double t916;
  double t918;
  double t922;
  double t928;
  double t932;
  double t940;
  double t941;
  double t964;
  double t976;
  double t979;
  double t984;
  double t985;
  double t986;
  double t992;
  double t997;
  double t998;
  double t1000;
  double t1004;
  double t1005;
  double t1009;
  double t1010;
  double t1012;
  double t1013;
  double t1014;
  double t1016;
  double t1017;
  double t1018;
  double t1019;
  double t1020;
  double t1041;
  double t1042;
  double t1043;
  double t1044;
  double t1045;
  double t999;
  double t1011;
  double t1015;
  double t1021;
  double t1022;
  double t1030;
  double t1031;
  double t1032;
  double t1033;
  double t1034;
  double t920;
  double t942;
  double t963;
  double t965;
  double t966;
  double t1053;
  double t1054;
  double t1055;
  double t1056;
  double t1057;
  double t1067;
  double t1068;
  double t1069;
  double t1026;
  double t1027;
  double t1028;
  double t863;
  double t883;
  double t885;
  double t1049;
  double t1050;
  double t1051;
  t831 = Cos(var1[5]);
  t810 = Cos(var1[6]);
  t816 = Sin(var1[5]);
  t832 = Sin(var1[6]);
  t844 = Cos(var1[2]);
  t830 = -1.*t810*t816;
  t834 = -1.*t831*t832;
  t840 = t830 + t834;
  t807 = Sin(var1[2]);
  t843 = t807*t840;
  t846 = t831*t810;
  t850 = -1.*t816*t832;
  t852 = t846 + t850;
  t854 = t844*t852;
  t861 = t843 + t854;
  t891 = -0.022659*t810;
  t893 = -0.007367999999999986*t832;
  t900 = t891 + t893;
  t905 = -1.*t810;
  t906 = 1. + t905;
  t909 = -0.16*t906;
  t910 = -0.167368*t810;
  t911 = 0.022659*t832;
  t914 = t909 + t910 + t911;
  t935 = -0.007367999999999986*t810;
  t937 = t935 + t911;
  t923 = 0.022659*t810;
  t925 = 0.007367999999999986*t832;
  t926 = t923 + t925;
  t943 = t810*t914;
  t949 = t900*t832;
  t869 = t844*t840;
  t874 = -1.*t831*t810;
  t877 = t816*t832;
  t878 = t874 + t877;
  t879 = t807*t878;
  t881 = t869 + t879;
  t952 = t810*t816;
  t953 = t831*t832;
  t954 = t952 + t953;
  t955 = t844*t954;
  t960 = t807*t852;
  t961 = t955 + t960;
  t977 = -1.*t807*t852;
  t978 = t869 + t977;
  t988 = -1.*t816*t900;
  t990 = t831*t914;
  t991 = t988 + t990;
  t994 = t831*t900;
  t995 = t816*t914;
  t996 = t994 + t995;
  t1006 = t831*t926;
  t1007 = -1.*t816*t937;
  t1008 = t1006 + t1007;
  t1001 = t816*t926;
  t1002 = t831*t937;
  t1003 = t1001 + t1002;
  t980 = -1.*t807*t840;
  t981 = t844*t878;
  t982 = t980 + t981;
  t947 = -1.*t810*t937;
  t950 = t926*t832;
  t951 = t943 + t947 + t949 + t950;
  t974 = -1.*t807*t954;
  t975 = t974 + t854;
  t902 = -1.*t810*t900;
  t916 = t914*t832;
  t918 = t902 + t916;
  t922 = t810*t900;
  t928 = t810*t926;
  t932 = -1.*t914*t832;
  t940 = t937*t832;
  t941 = t922 + t928 + t932 + t940;
  t964 = t943 + t949;
  t976 = 0.14994*t861*t975;
  t979 = 0.14994*t978*t961;
  t984 = 0.14994*t861*t982;
  t985 = 0.14994*t978*t881;
  t986 = t976 + t979 + t984 + t985;
  t992 = -1.*t991*t840;
  t997 = -1.*t996*t852;
  t998 = t992 + t997;
  t1000 = t991*t840;
  t1004 = t1003*t954;
  t1005 = t996*t852;
  t1009 = t1008*t852;
  t1010 = t1000 + t1004 + t1005 + t1009;
  t1012 = t996*t954;
  t1013 = t991*t852;
  t1014 = t1012 + t1013;
  t1016 = -1.*t996*t840;
  t1017 = -1.*t1008*t840;
  t1018 = -1.*t1003*t852;
  t1019 = -1.*t991*t878;
  t1020 = t1016 + t1017 + t1018 + t1019;
  t1041 = 0.14994*t978*t998;
  t1042 = 0.14994*t978*t1010;
  t1043 = 0.14994*t1014*t982;
  t1044 = 0.14994*t975*t1020;
  t1045 = t1041 + t1042 + t1043 + t1044;
  t999 = 0.14994*t861*t998;
  t1011 = 0.14994*t861*t1010;
  t1015 = 0.14994*t1014*t881;
  t1021 = 0.14994*t961*t1020;
  t1022 = t999 + t1011 + t1015 + t1021;
  t1030 = 0.14994*t951*t975;
  t1031 = 0.14994*t918*t978;
  t1032 = 0.14994*t941*t978;
  t1033 = 0.14994*t964*t982;
  t1034 = t1030 + t1031 + t1032 + t1033;
  t920 = 0.14994*t918*t861;
  t942 = 0.14994*t941*t861;
  t963 = 0.14994*t951*t961;
  t965 = 0.14994*t964*t881;
  t966 = t920 + t942 + t963 + t965;
  t1053 = 0.14994*t941*t1014;
  t1054 = 0.14994*t951*t998;
  t1055 = 0.14994*t964*t1010;
  t1056 = 0.14994*t918*t1020;
  t1057 = t1053 + t1054 + t1055 + t1056;
  t1067 = 0.0033974904599999994*t951;
  t1068 = -0.0011047579199999977*t941;
  t1069 = t1067 + t1068;
  t1026 = 0.0033974904599999994*t978;
  t1027 = -0.0011047579199999977*t982;
  t1028 = t1026 + t1027;
  t863 = 0.0033974904599999994*t861;
  t883 = -0.0011047579199999977*t881;
  t885 = t863 + t883;
  t1049 = -0.0011047579199999977*t1010;
  t1050 = 0.0033974904599999994*t1020;
  t1051 = t1049 + t1050;
  p_output1[0]=var2[6]*(-0.5*(0.29988*t861*t881 + 0.29988*t861*t961)*var2[0] - 0.5*t986*var2[1] - 0.5*t1022*var2[2] - 0.5*t966*var2[5] - 0.5*t885*var2[6]);
  p_output1[1]=var2[6]*(-0.5*t986*var2[0] - 0.5*(0.29988*t975*t978 + 0.29988*t978*t982)*var2[1] - 0.5*t1045*var2[2] - 0.5*t1034*var2[5] - 0.5*t1028*var2[6]);
  p_output1[2]=var2[6]*(-0.5*t1022*var2[0] - 0.5*t1045*var2[1] - 0.5*(0.29988*t1010*t1014 + 0.29988*t1020*t998)*var2[2] - 0.5*t1057*var2[5] - 0.5*t1051*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[6]*(-0.5*t966*var2[0] - 0.5*t1034*var2[1] - 0.5*t1057*var2[2] - 0.5*(0.29988*t918*t951 + 0.29988*t941*t964)*var2[5] - 0.5*t1069*var2[6]);
  p_output1[6]=(-0.5*t885*var2[0] - 0.5*t1028*var2[1] - 0.5*t1051*var2[2] - 0.5*t1069*var2[5])*var2[6];
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

#include "Ce2_vec_L5_J7_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L5_J7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
