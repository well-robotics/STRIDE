/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:32:27 GMT-05:00
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
  double t969;
  double t943;
  double t947;
  double t970;
  double t980;
  double t957;
  double t972;
  double t974;
  double t939;
  double t1025;
  double t1026;
  double t977;
  double t985;
  double t990;
  double t994;
  double t999;
  double t1000;
  double t1028;
  double t1033;
  double t1034;
  double t1035;
  double t1042;
  double t1044;
  double t1045;
  double t1046;
  double t1070;
  double t1073;
  double t1060;
  double t1061;
  double t1063;
  double t1079;
  double t1082;
  double t1084;
  double t1085;
  double t1006;
  double t1008;
  double t1009;
  double t1010;
  double t1015;
  double t1017;
  double t1089;
  double t1093;
  double t1094;
  double t1095;
  double t1096;
  double t1097;
  double t1116;
  double t1118;
  double t1127;
  double t1129;
  double t1130;
  double t1131;
  double t1132;
  double t1133;
  double t1135;
  double t1137;
  double t1139;
  double t1140;
  double t1141;
  double t1142;
  double t1152;
  double t1153;
  double t1154;
  double t1147;
  double t1148;
  double t1149;
  double t1120;
  double t1121;
  double t1122;
  double t1083;
  double t1087;
  double t1088;
  double t1112;
  double t1114;
  double t1021;
  double t1039;
  double t1040;
  double t1050;
  double t1053;
  double t1056;
  double t1057;
  double t1064;
  double t1066;
  double t1069;
  double t1076;
  double t1077;
  double t1099;
  double t1115;
  double t1119;
  double t1123;
  double t1124;
  double t1125;
  double t1136;
  double t1143;
  double t1144;
  double t1146;
  double t1150;
  double t1151;
  double t1155;
  double t1156;
  double t1158;
  double t1159;
  double t1160;
  double t1162;
  double t1163;
  double t1164;
  double t1165;
  double t1166;
  double t1187;
  double t1188;
  double t1189;
  double t1190;
  double t1191;
  double t1145;
  double t1157;
  double t1161;
  double t1167;
  double t1168;
  double t1176;
  double t1177;
  double t1178;
  double t1179;
  double t1180;
  double t1054;
  double t1078;
  double t1098;
  double t1100;
  double t1105;
  double t1199;
  double t1200;
  double t1201;
  double t1202;
  double t1203;
  double t1213;
  double t1214;
  double t1215;
  double t1172;
  double t1173;
  double t1174;
  double t1002;
  double t1018;
  double t1019;
  double t1195;
  double t1196;
  double t1197;
  t969 = Cos(var1[5]);
  t943 = Cos(var1[6]);
  t947 = Sin(var1[5]);
  t970 = Sin(var1[6]);
  t980 = Cos(var1[2]);
  t957 = -1.*t943*t947;
  t972 = -1.*t969*t970;
  t974 = t957 + t972;
  t939 = Sin(var1[2]);
  t1025 = -1.*t943;
  t1026 = 1. + t1025;
  t977 = t939*t974;
  t985 = t969*t943;
  t990 = -1.*t947*t970;
  t994 = t985 + t990;
  t999 = t980*t994;
  t1000 = t977 + t999;
  t1028 = -0.0265*t1026;
  t1033 = -0.025226*t943;
  t1034 = -0.07700600000000002*t970;
  t1035 = t1028 + t1033 + t1034;
  t1042 = -0.2375*t1026;
  t1044 = -0.314506*t943;
  t1045 = -0.0012740000000000008*t970;
  t1046 = t1042 + t1044 + t1045;
  t1070 = -0.07700600000000002*t943;
  t1073 = t1070 + t1045;
  t1060 = -0.0012740000000000008*t943;
  t1061 = 0.07700600000000002*t970;
  t1063 = t1060 + t1061;
  t1079 = 0.0695*t943;
  t1082 = t943*t1046;
  t1084 = 0.0265*t970;
  t1085 = t1035*t970;
  t1006 = t980*t974;
  t1008 = -1.*t969*t943;
  t1009 = t947*t970;
  t1010 = t1008 + t1009;
  t1015 = t939*t1010;
  t1017 = t1006 + t1015;
  t1089 = t943*t947;
  t1093 = t969*t970;
  t1094 = t1089 + t1093;
  t1095 = t980*t1094;
  t1096 = t939*t994;
  t1097 = t1095 + t1096;
  t1116 = -1.*t939*t994;
  t1118 = t1006 + t1116;
  t1127 = -1.*t969;
  t1129 = 1. + t1127;
  t1130 = -0.0695*t1129;
  t1131 = -0.0265*t947;
  t1132 = -1.*t947*t1035;
  t1133 = t969*t1046;
  t1135 = t1130 + t1131 + t1132 + t1133;
  t1137 = -0.0265*t1129;
  t1139 = 0.0695*t947;
  t1140 = t969*t1035;
  t1141 = t947*t1046;
  t1142 = t1137 + t1139 + t1140 + t1141;
  t1152 = -1.*t947*t1073;
  t1153 = t969*t1063;
  t1154 = t1152 + t1153;
  t1147 = t969*t1073;
  t1148 = t947*t1063;
  t1149 = t1147 + t1148;
  t1120 = -1.*t939*t974;
  t1121 = t980*t1010;
  t1122 = t1120 + t1121;
  t1083 = -1.*t943*t1073;
  t1087 = t1063*t970;
  t1088 = t1079 + t1082 + t1083 + t1084 + t1085 + t1087;
  t1112 = -1.*t939*t1094;
  t1114 = t1112 + t999;
  t1021 = -0.0265*t943;
  t1039 = -1.*t943*t1035;
  t1040 = 0.0695*t970;
  t1050 = t1046*t970;
  t1053 = t1021 + t1039 + t1040 + t1050;
  t1056 = 0.0265*t943;
  t1057 = t943*t1035;
  t1064 = t943*t1063;
  t1066 = -0.0695*t970;
  t1069 = -1.*t1046*t970;
  t1076 = t1073*t970;
  t1077 = t1056 + t1057 + t1064 + t1066 + t1069 + t1076;
  t1099 = t1079 + t1082 + t1084 + t1085;
  t1115 = 0.19964*t1000*t1114;
  t1119 = 0.19964*t1118*t1097;
  t1123 = 0.19964*t1000*t1122;
  t1124 = 0.19964*t1118*t1017;
  t1125 = t1115 + t1119 + t1123 + t1124;
  t1136 = -1.*t1135*t974;
  t1143 = -1.*t1142*t994;
  t1144 = t1136 + t1143;
  t1146 = t1135*t974;
  t1150 = t1149*t1094;
  t1151 = t1142*t994;
  t1155 = t1154*t994;
  t1156 = t1146 + t1150 + t1151 + t1155;
  t1158 = t1142*t1094;
  t1159 = t1135*t994;
  t1160 = t1158 + t1159;
  t1162 = -1.*t1142*t974;
  t1163 = -1.*t1154*t974;
  t1164 = -1.*t1149*t994;
  t1165 = -1.*t1135*t1010;
  t1166 = t1162 + t1163 + t1164 + t1165;
  t1187 = 0.19964*t1118*t1144;
  t1188 = 0.19964*t1118*t1156;
  t1189 = 0.19964*t1160*t1122;
  t1190 = 0.19964*t1114*t1166;
  t1191 = t1187 + t1188 + t1189 + t1190;
  t1145 = 0.19964*t1000*t1144;
  t1157 = 0.19964*t1000*t1156;
  t1161 = 0.19964*t1160*t1017;
  t1167 = 0.19964*t1097*t1166;
  t1168 = t1145 + t1157 + t1161 + t1167;
  t1176 = 0.19964*t1088*t1114;
  t1177 = 0.19964*t1053*t1118;
  t1178 = 0.19964*t1077*t1118;
  t1179 = 0.19964*t1099*t1122;
  t1180 = t1176 + t1177 + t1178 + t1179;
  t1054 = 0.19964*t1053*t1000;
  t1078 = 0.19964*t1077*t1000;
  t1098 = 0.19964*t1088*t1097;
  t1100 = 0.19964*t1099*t1017;
  t1105 = t1054 + t1078 + t1098 + t1100;
  t1199 = 0.19964*t1077*t1160;
  t1200 = 0.19964*t1088*t1144;
  t1201 = 0.19964*t1099*t1156;
  t1202 = 0.19964*t1053*t1166;
  t1203 = t1199 + t1200 + t1201 + t1202;
  t1213 = -0.015373477840000005*t1077;
  t1214 = -0.0002543413600000002*t1088;
  t1215 = t1213 + t1214;
  t1172 = -0.0002543413600000002*t1118;
  t1173 = -0.015373477840000005*t1122;
  t1174 = t1172 + t1173;
  t1002 = -0.0002543413600000002*t1000;
  t1018 = -0.015373477840000005*t1017;
  t1019 = t1002 + t1018;
  t1195 = -0.015373477840000005*t1156;
  t1196 = -0.0002543413600000002*t1166;
  t1197 = t1195 + t1196;
  p_output1[0]=var2[6]*(-0.5*(0.39928*t1000*t1017 + 0.39928*t1000*t1097)*var2[0] - 0.5*t1125*var2[1] - 0.5*t1168*var2[2] - 0.5*t1105*var2[5] - 0.5*t1019*var2[6]);
  p_output1[1]=var2[6]*(-0.5*t1125*var2[0] - 0.5*(0.39928*t1114*t1118 + 0.39928*t1118*t1122)*var2[1] - 0.5*t1191*var2[2] - 0.5*t1180*var2[5] - 0.5*t1174*var2[6]);
  p_output1[2]=var2[6]*(-0.5*t1168*var2[0] - 0.5*t1191*var2[1] - 0.5*(0.39928*t1156*t1160 + 0.39928*t1144*t1166)*var2[2] - 0.5*t1203*var2[5] - 0.5*t1197*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[6]*(-0.5*t1105*var2[0] - 0.5*t1180*var2[1] - 0.5*t1203*var2[2] - 0.5*(0.39928*t1053*t1088 + 0.39928*t1077*t1099)*var2[5] - 0.5*t1215*var2[6]);
  p_output1[6]=(-0.5*t1019*var2[0] - 0.5*t1174*var2[1] - 0.5*t1197*var2[2] - 0.5*t1215*var2[5])*var2[6];
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

#include "Ce2_vec_L5_J7_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L5_J7_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
