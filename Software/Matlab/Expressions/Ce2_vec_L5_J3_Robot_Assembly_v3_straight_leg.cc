/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:32:17 GMT-05:00
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
  double t882;
  double t869;
  double t874;
  double t886;
  double t892;
  double t868;
  double t897;
  double t899;
  double t902;
  double t948;
  double t949;
  double t881;
  double t887;
  double t888;
  double t890;
  double t909;
  double t918;
  double t970;
  double t972;
  double t973;
  double t974;
  double t953;
  double t956;
  double t957;
  double t959;
  double t921;
  double t924;
  double t926;
  double t927;
  double t928;
  double t934;
  double t1000;
  double t1002;
  double t1005;
  double t1006;
  double t1008;
  double t1029;
  double t1031;
  double t1039;
  double t1040;
  double t1041;
  double t1042;
  double t1043;
  double t1032;
  double t1033;
  double t1034;
  double t1035;
  double t1037;
  double t1012;
  double t1014;
  double t1015;
  double t1023;
  double t1025;
  double t985;
  double t989;
  double t990;
  double t994;
  double t995;
  double t947;
  double t964;
  double t969;
  double t976;
  double t977;
  double t1017;
  double t1018;
  double t1019;
  double t1020;
  double t1021;
  double t1026;
  double t1027;
  double t1038;
  double t1044;
  double t1045;
  double t1047;
  double t1049;
  double t1050;
  double t1046;
  double t1051;
  double t1052;
  double t1069;
  double t1070;
  double t1071;
  double t980;
  double t996;
  double t998;
  double t1060;
  double t1061;
  double t1062;
  double t920;
  double t939;
  double t941;
  double t1056;
  double t1057;
  double t1058;
  t882 = Cos(var1[5]);
  t869 = Cos(var1[6]);
  t874 = Sin(var1[5]);
  t886 = Sin(var1[6]);
  t892 = Cos(var1[2]);
  t868 = Sin(var1[2]);
  t897 = t882*t869;
  t899 = -1.*t874*t886;
  t902 = t897 + t899;
  t948 = -1.*t869;
  t949 = 1. + t948;
  t881 = t869*t874;
  t887 = t882*t886;
  t888 = t881 + t887;
  t890 = -1.*t868*t888;
  t909 = t892*t902;
  t918 = t890 + t909;
  t970 = -0.2375*t949;
  t972 = -0.314506*t869;
  t973 = -0.0012740000000000008*t886;
  t974 = t970 + t972 + t973;
  t953 = -0.0265*t949;
  t956 = -0.025226*t869;
  t957 = -0.07700600000000002*t886;
  t959 = t953 + t956 + t957;
  t921 = -1.*t869*t874;
  t924 = -1.*t882*t886;
  t926 = t921 + t924;
  t927 = t892*t926;
  t928 = -1.*t868*t902;
  t934 = t927 + t928;
  t1000 = t868*t926;
  t1002 = t1000 + t909;
  t1005 = t892*t888;
  t1006 = t868*t902;
  t1008 = t1005 + t1006;
  t1029 = -1.*t882;
  t1031 = 1. + t1029;
  t1039 = -0.0695*t1031;
  t1040 = -0.0265*t874;
  t1041 = -1.*t874*t959;
  t1042 = t882*t974;
  t1043 = t1039 + t1040 + t1041 + t1042;
  t1032 = -0.0265*t1031;
  t1033 = 0.0695*t874;
  t1034 = t882*t959;
  t1035 = t874*t974;
  t1037 = t1032 + t1033 + t1034 + t1035;
  t1012 = -1.*t868*t926;
  t1014 = -1.*t892*t902;
  t1015 = t1012 + t1014;
  t1023 = -1.*t892*t888;
  t1025 = t1023 + t928;
  t985 = 0.0695*t869;
  t989 = t869*t974;
  t990 = 0.0265*t886;
  t994 = t959*t886;
  t995 = t985 + t989 + t990 + t994;
  t947 = -0.0265*t869;
  t964 = -1.*t869*t959;
  t969 = 0.0695*t886;
  t976 = t974*t886;
  t977 = t947 + t964 + t969 + t976;
  t1017 = 0.19964*t1015*t1002;
  t1018 = Power(t918,2);
  t1019 = 0.19964*t1018;
  t1020 = Power(t934,2);
  t1021 = 0.19964*t1020;
  t1026 = 0.19964*t1025*t1008;
  t1027 = t1017 + t1019 + t1021 + t1026;
  t1038 = t1037*t888;
  t1044 = t1043*t902;
  t1045 = t1038 + t1044;
  t1047 = -1.*t1043*t926;
  t1049 = -1.*t1037*t902;
  t1050 = t1047 + t1049;
  t1046 = 0.19964*t934*t1045;
  t1051 = 0.19964*t918*t1050;
  t1052 = t1046 + t1051;
  t1069 = 0.19964*t1015*t1045;
  t1070 = 0.19964*t1025*t1050;
  t1071 = t1069 + t1070;
  t980 = 0.19964*t977*t918;
  t996 = 0.19964*t995*t934;
  t998 = t980 + t996;
  t1060 = 0.19964*t995*t1015;
  t1061 = 0.19964*t977*t1025;
  t1062 = t1060 + t1061;
  t920 = -0.0002543413600000002*t918;
  t939 = -0.015373477840000005*t934;
  t941 = t920 + t939;
  t1056 = -0.015373477840000005*t1015;
  t1057 = -0.0002543413600000002*t1025;
  t1058 = t1056 + t1057;
  p_output1[0]=var2[2]*(-0.5*(0.39928*t1008*t918 + 0.39928*t1002*t934)*var2[0] - 0.5*t1027*var2[1] - 0.5*t1052*var2[2] - 0.5*t998*var2[5] - 0.5*t941*var2[6]);
  p_output1[1]=var2[2]*(-0.5*t1027*var2[0] - 0.5*(0.39928*t1025*t918 + 0.39928*t1015*t934)*var2[1] - 0.5*t1071*var2[2] - 0.5*t1062*var2[5] - 0.5*t1058*var2[6]);
  p_output1[2]=(-0.5*t1052*var2[0] - 0.5*t1071*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t998*var2[0] - 0.5*t1062*var2[1])*var2[2];
  p_output1[6]=(-0.5*t941*var2[0] - 0.5*t1058*var2[1])*var2[2];
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

#include "Ce2_vec_L5_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L5_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
