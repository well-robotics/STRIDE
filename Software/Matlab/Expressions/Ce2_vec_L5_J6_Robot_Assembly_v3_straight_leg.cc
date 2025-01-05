/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:32:25 GMT-05:00
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
  double t918;
  double t887;
  double t897;
  double t920;
  double t943;
  double t909;
  double t927;
  double t934;
  double t881;
  double t995;
  double t996;
  double t939;
  double t947;
  double t953;
  double t956;
  double t957;
  double t969;
  double t1008;
  double t1009;
  double t1010;
  double t1011;
  double t999;
  double t1000;
  double t1002;
  double t1003;
  double t972;
  double t973;
  double t974;
  double t976;
  double t977;
  double t980;
  double t1034;
  double t1035;
  double t1038;
  double t1039;
  double t1040;
  double t1042;
  double t1055;
  double t1056;
  double t1067;
  double t1068;
  double t1069;
  double t1070;
  double t1072;
  double t1073;
  double t1074;
  double t1076;
  double t1077;
  double t1078;
  double t1079;
  double t1080;
  double t1088;
  double t1089;
  double t1090;
  double t1091;
  double t1092;
  double t1085;
  double t1086;
  double t1059;
  double t1060;
  double t1061;
  double t994;
  double t1005;
  double t1006;
  double t1015;
  double t1017;
  double t1019;
  double t1020;
  double t1021;
  double t1025;
  double t1026;
  double t1051;
  double t1053;
  double t1054;
  double t1057;
  double t1063;
  double t1064;
  double t1065;
  double t1075;
  double t1081;
  double t1082;
  double t1084;
  double t1087;
  double t1093;
  double t1094;
  double t1095;
  double t1097;
  double t1098;
  double t1099;
  double t1101;
  double t1102;
  double t1103;
  double t1104;
  double t1105;
  double t1124;
  double t1125;
  double t1126;
  double t1127;
  double t1128;
  double t1083;
  double t1096;
  double t1100;
  double t1106;
  double t1107;
  double t1115;
  double t1116;
  double t1117;
  double t1018;
  double t1028;
  double t1032;
  double t1136;
  double t1137;
  double t1138;
  double t1111;
  double t1112;
  double t1113;
  double t970;
  double t985;
  double t989;
  double t1132;
  double t1133;
  double t1134;
  t918 = Cos(var1[5]);
  t887 = Cos(var1[6]);
  t897 = Sin(var1[5]);
  t920 = Sin(var1[6]);
  t943 = Cos(var1[2]);
  t909 = -1.*t887*t897;
  t927 = -1.*t918*t920;
  t934 = t909 + t927;
  t881 = Sin(var1[2]);
  t995 = -1.*t887;
  t996 = 1. + t995;
  t939 = t881*t934;
  t947 = t918*t887;
  t953 = -1.*t897*t920;
  t956 = t947 + t953;
  t957 = t943*t956;
  t969 = t939 + t957;
  t1008 = -0.2375*t996;
  t1009 = -0.314506*t887;
  t1010 = -0.0012740000000000008*t920;
  t1011 = t1008 + t1009 + t1010;
  t999 = -0.0265*t996;
  t1000 = -0.025226*t887;
  t1002 = -0.07700600000000002*t920;
  t1003 = t999 + t1000 + t1002;
  t972 = t943*t934;
  t973 = -1.*t918*t887;
  t974 = t897*t920;
  t976 = t973 + t974;
  t977 = t881*t976;
  t980 = t972 + t977;
  t1034 = t887*t897;
  t1035 = t918*t920;
  t1038 = t1034 + t1035;
  t1039 = t943*t1038;
  t1040 = t881*t956;
  t1042 = t1039 + t1040;
  t1055 = -1.*t881*t956;
  t1056 = t972 + t1055;
  t1067 = -1.*t918;
  t1068 = 1. + t1067;
  t1069 = -0.0695*t1068;
  t1070 = -0.0265*t897;
  t1072 = -1.*t897*t1003;
  t1073 = t918*t1011;
  t1074 = t1069 + t1070 + t1072 + t1073;
  t1076 = -0.0265*t1068;
  t1077 = 0.0695*t897;
  t1078 = t918*t1003;
  t1079 = t897*t1011;
  t1080 = t1076 + t1077 + t1078 + t1079;
  t1088 = -0.0265*t918;
  t1089 = -0.0695*t897;
  t1090 = -1.*t918*t1003;
  t1091 = -1.*t897*t1011;
  t1092 = t1088 + t1089 + t1090 + t1091;
  t1085 = 0.0695*t918;
  t1086 = t1085 + t1070 + t1072 + t1073;
  t1059 = -1.*t881*t934;
  t1060 = t943*t976;
  t1061 = t1059 + t1060;
  t994 = -0.0265*t887;
  t1005 = -1.*t887*t1003;
  t1006 = 0.0695*t920;
  t1015 = t1011*t920;
  t1017 = t994 + t1005 + t1006 + t1015;
  t1019 = 0.0695*t887;
  t1020 = t887*t1011;
  t1021 = 0.0265*t920;
  t1025 = t1003*t920;
  t1026 = t1019 + t1020 + t1021 + t1025;
  t1051 = -1.*t881*t1038;
  t1053 = t1051 + t957;
  t1054 = 0.19964*t969*t1053;
  t1057 = 0.19964*t1056*t1042;
  t1063 = 0.19964*t969*t1061;
  t1064 = 0.19964*t1056*t980;
  t1065 = t1054 + t1057 + t1063 + t1064;
  t1075 = -1.*t1074*t934;
  t1081 = -1.*t1080*t956;
  t1082 = t1075 + t1081;
  t1084 = t1074*t934;
  t1087 = t1086*t1038;
  t1093 = t1092*t956;
  t1094 = t1080*t956;
  t1095 = t1084 + t1087 + t1093 + t1094;
  t1097 = t1080*t1038;
  t1098 = t1074*t956;
  t1099 = t1097 + t1098;
  t1101 = -1.*t1092*t934;
  t1102 = -1.*t1080*t934;
  t1103 = -1.*t1086*t956;
  t1104 = -1.*t1074*t976;
  t1105 = t1101 + t1102 + t1103 + t1104;
  t1124 = 0.19964*t1056*t1082;
  t1125 = 0.19964*t1056*t1095;
  t1126 = 0.19964*t1099*t1061;
  t1127 = 0.19964*t1053*t1105;
  t1128 = t1124 + t1125 + t1126 + t1127;
  t1083 = 0.19964*t969*t1082;
  t1096 = 0.19964*t969*t1095;
  t1100 = 0.19964*t1099*t980;
  t1106 = 0.19964*t1042*t1105;
  t1107 = t1083 + t1096 + t1100 + t1106;
  t1115 = 0.19964*t1017*t1056;
  t1116 = 0.19964*t1026*t1061;
  t1117 = t1115 + t1116;
  t1018 = 0.19964*t1017*t969;
  t1028 = 0.19964*t1026*t980;
  t1032 = t1018 + t1028;
  t1136 = 0.19964*t1026*t1095;
  t1137 = 0.19964*t1017*t1105;
  t1138 = t1136 + t1137;
  t1111 = -0.0002543413600000002*t1056;
  t1112 = -0.015373477840000005*t1061;
  t1113 = t1111 + t1112;
  t970 = -0.0002543413600000002*t969;
  t985 = -0.015373477840000005*t980;
  t989 = t970 + t985;
  t1132 = -0.015373477840000005*t1095;
  t1133 = -0.0002543413600000002*t1105;
  t1134 = t1132 + t1133;
  p_output1[0]=var2[5]*(-0.5*(0.39928*t1042*t969 + 0.39928*t969*t980)*var2[0] - 0.5*t1065*var2[1] - 0.5*t1107*var2[2] - 0.5*t1032*var2[5] - 0.5*t989*var2[6]);
  p_output1[1]=var2[5]*(-0.5*t1065*var2[0] - 0.5*(0.39928*t1053*t1056 + 0.39928*t1056*t1061)*var2[1] - 0.5*t1128*var2[2] - 0.5*t1117*var2[5] - 0.5*t1113*var2[6]);
  p_output1[2]=var2[5]*(-0.5*t1107*var2[0] - 0.5*t1128*var2[1] - 0.5*(0.39928*t1095*t1099 + 0.39928*t1082*t1105)*var2[2] - 0.5*t1138*var2[5] - 0.5*t1134*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t1032*var2[0] - 0.5*t1117*var2[1] - 0.5*t1138*var2[2])*var2[5];
  p_output1[6]=(-0.5*t989*var2[0] - 0.5*t1113*var2[1] - 0.5*t1134*var2[2])*var2[5];
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

#include "Ce2_vec_L5_J6_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L5_J6_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
