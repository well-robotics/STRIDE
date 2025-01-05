/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:31:48 GMT-05:00
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
  double t812;
  double t787;
  double t797;
  double t814;
  double t828;
  double t808;
  double t816;
  double t821;
  double t731;
  double t869;
  double t861;
  double t862;
  double t830;
  double t834;
  double t835;
  double t888;
  double t890;
  double t892;
  double t893;
  double t864;
  double t868;
  double t870;
  double t825;
  double t840;
  double t841;
  double t882;
  double t883;
  double t886;
  double t876;
  double t879;
  double t859;
  double t874;
  double t881;
  double t896;
  double t845;
  double t846;
  double t847;
  double t849;
  double t851;
  double t852;
  double t899;
  double t902;
  double t903;
  double t906;
  double t909;
  double t912;
  double t945;
  double t947;
  double t960;
  double t962;
  double t963;
  double t964;
  double t965;
  double t966;
  double t968;
  double t970;
  double t972;
  double t973;
  double t974;
  double t975;
  double t985;
  double t986;
  double t987;
  double t980;
  double t981;
  double t982;
  double t955;
  double t956;
  double t918;
  double t919;
  double t920;
  double t921;
  double t924;
  double t926;
  double t927;
  double t928;
  double t929;
  double t932;
  double t933;
  double t934;
  double t880;
  double t887;
  double t897;
  double t949;
  double t951;
  double t952;
  double t936;
  double t948;
  double t953;
  double t954;
  double t957;
  double t958;
  double t969;
  double t976;
  double t977;
  double t979;
  double t983;
  double t984;
  double t988;
  double t989;
  double t991;
  double t992;
  double t993;
  double t995;
  double t996;
  double t997;
  double t998;
  double t999;
  double t978;
  double t990;
  double t994;
  double t1000;
  double t1001;
  double t1020;
  double t1021;
  double t1022;
  double t1023;
  double t1024;
  double t913;
  double t925;
  double t935;
  double t937;
  double t938;
  double t1009;
  double t1010;
  double t1011;
  double t1012;
  double t1013;
  double t1032;
  double t1033;
  double t1034;
  double t1035;
  double t1036;
  double t1046;
  double t1047;
  double t1048;
  double t842;
  double t856;
  double t857;
  double t1005;
  double t1006;
  double t1007;
  double t1028;
  double t1029;
  double t1030;
  t812 = Cos(var1[3]);
  t787 = Cos(var1[4]);
  t797 = Sin(var1[3]);
  t814 = Sin(var1[4]);
  t828 = Sin(var1[2]);
  t808 = -1.*t787*t797;
  t816 = -1.*t812*t814;
  t821 = t808 + t816;
  t731 = Cos(var1[2]);
  t869 = 0.07701400000000003*t814;
  t861 = -1.*t787;
  t862 = 1. + t861;
  t830 = t812*t787;
  t834 = -1.*t797*t814;
  t835 = t830 + t834;
  t888 = -0.2375*t862;
  t890 = -0.314514*t787;
  t892 = 0.0012709999999999978*t814;
  t893 = t888 + t890 + t892;
  t864 = -0.0265*t862;
  t868 = -0.025229*t787;
  t870 = t864 + t868 + t869;
  t825 = t731*t821;
  t840 = t828*t835;
  t841 = t825 + t840;
  t882 = 0.07701400000000003*t787;
  t883 = -0.0012709999999999978*t814;
  t886 = t882 + t883;
  t876 = 0.0012709999999999978*t787;
  t879 = t876 + t869;
  t859 = 0.0265*t787;
  t874 = t787*t870;
  t881 = 0.0695*t814;
  t896 = t893*t814;
  t845 = t828*t821;
  t846 = -1.*t812*t787;
  t847 = t797*t814;
  t849 = t846 + t847;
  t851 = t731*t849;
  t852 = t845 + t851;
  t899 = t787*t797;
  t902 = t812*t814;
  t903 = t899 + t902;
  t906 = t828*t903;
  t909 = t731*t835;
  t912 = t906 + t909;
  t945 = -1.*t828*t821;
  t947 = t945 + t909;
  t960 = -1.*t812;
  t962 = 1. + t960;
  t963 = -0.0265*t962;
  t964 = -0.0695*t797;
  t965 = -1.*t797*t893;
  t966 = t812*t870;
  t968 = t963 + t964 + t965 + t966;
  t970 = -0.0695*t962;
  t972 = 0.0265*t797;
  t973 = t812*t893;
  t974 = t797*t870;
  t975 = t970 + t972 + t973 + t974;
  t985 = t812*t886;
  t986 = -1.*t797*t879;
  t987 = t985 + t986;
  t980 = t797*t886;
  t981 = t812*t879;
  t982 = t980 + t981;
  t955 = -1.*t828*t849;
  t956 = t825 + t955;
  t918 = -0.0695*t787;
  t919 = -1.*t787*t893;
  t920 = 0.0265*t814;
  t921 = t870*t814;
  t924 = t918 + t919 + t920 + t921;
  t926 = 0.0695*t787;
  t927 = t787*t886;
  t928 = t787*t893;
  t929 = -0.0265*t814;
  t932 = -1.*t870*t814;
  t933 = t879*t814;
  t934 = t926 + t927 + t928 + t929 + t932 + t933;
  t880 = -1.*t787*t879;
  t887 = t886*t814;
  t897 = t859 + t874 + t880 + t881 + t887 + t896;
  t949 = t731*t903;
  t951 = -1.*t828*t835;
  t952 = t949 + t951;
  t936 = t859 + t874 + t881 + t896;
  t948 = 0.19964*t947*t912;
  t953 = 0.19964*t952*t841;
  t954 = 0.19964*t947*t852;
  t957 = 0.19964*t841*t956;
  t958 = t948 + t953 + t954 + t957;
  t969 = t968*t821;
  t976 = t975*t835;
  t977 = t969 + t976;
  t979 = -1.*t968*t821;
  t983 = -1.*t982*t903;
  t984 = -1.*t975*t835;
  t988 = -1.*t987*t835;
  t989 = t979 + t983 + t984 + t988;
  t991 = -1.*t975*t903;
  t992 = -1.*t968*t835;
  t993 = t991 + t992;
  t995 = t975*t821;
  t996 = t987*t821;
  t997 = t982*t835;
  t998 = t968*t849;
  t999 = t995 + t996 + t997 + t998;
  t978 = 0.19964*t841*t977;
  t990 = 0.19964*t841*t989;
  t994 = 0.19964*t993*t852;
  t1000 = 0.19964*t912*t999;
  t1001 = t978 + t990 + t994 + t1000;
  t1020 = 0.19964*t947*t977;
  t1021 = 0.19964*t947*t989;
  t1022 = 0.19964*t993*t956;
  t1023 = 0.19964*t952*t999;
  t1024 = t1020 + t1021 + t1022 + t1023;
  t913 = 0.19964*t897*t912;
  t925 = 0.19964*t924*t841;
  t935 = 0.19964*t934*t841;
  t937 = 0.19964*t936*t852;
  t938 = t913 + t925 + t935 + t937;
  t1009 = 0.19964*t924*t947;
  t1010 = 0.19964*t934*t947;
  t1011 = 0.19964*t897*t952;
  t1012 = 0.19964*t936*t956;
  t1013 = t1009 + t1010 + t1011 + t1012;
  t1032 = 0.19964*t934*t993;
  t1033 = 0.19964*t897*t977;
  t1034 = 0.19964*t936*t989;
  t1035 = 0.19964*t924*t999;
  t1036 = t1032 + t1033 + t1034 + t1035;
  t1046 = 0.015375074960000006*t897;
  t1047 = 0.0002537424399999996*t934;
  t1048 = t1046 + t1047;
  t842 = 0.015375074960000006*t841;
  t856 = 0.0002537424399999996*t852;
  t857 = t842 + t856;
  t1005 = 0.015375074960000006*t947;
  t1006 = 0.0002537424399999996*t956;
  t1007 = t1005 + t1006;
  t1028 = 0.0002537424399999996*t989;
  t1029 = 0.015375074960000006*t999;
  t1030 = t1028 + t1029;
  p_output1[0]=var2[4]*(-0.5*(0.39928*t841*t852 + 0.39928*t841*t912)*var2[0] - 0.5*t958*var2[1] - 0.5*t1001*var2[2] - 0.5*t938*var2[3] - 0.5*t857*var2[4]);
  p_output1[1]=var2[4]*(-0.5*t958*var2[0] - 0.5*(0.39928*t947*t952 + 0.39928*t947*t956)*var2[1] - 0.5*t1024*var2[2] - 0.5*t1013*var2[3] - 0.5*t1007*var2[4]);
  p_output1[2]=var2[4]*(-0.5*t1001*var2[0] - 0.5*t1024*var2[1] - 0.5*(0.39928*t989*t993 + 0.39928*t977*t999)*var2[2] - 0.5*t1036*var2[3] - 0.5*t1030*var2[4]);
  p_output1[3]=var2[4]*(-0.5*t938*var2[0] - 0.5*t1013*var2[1] - 0.5*t1036*var2[2] - 0.5*(0.39928*t897*t924 + 0.39928*t934*t936)*var2[3] - 0.5*t1048*var2[4]);
  p_output1[4]=(-0.5*t857*var2[0] - 0.5*t1007*var2[1] - 0.5*t1030*var2[2] - 0.5*t1048*var2[3])*var2[4];
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

#include "Ce2_vec_L3_J5_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L3_J5_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
