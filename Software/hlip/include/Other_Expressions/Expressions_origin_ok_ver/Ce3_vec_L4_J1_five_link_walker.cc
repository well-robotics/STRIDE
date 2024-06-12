/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:49 GMT-05:00
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
  double t1068;
  double t1053;
  double t1058;
  double t1074;
  double t1087;
  double t1050;
  double t1088;
  double t1089;
  double t1090;
  double t1066;
  double t1081;
  double t1083;
  double t1084;
  double t1091;
  double t1092;
  double t1107;
  double t1108;
  double t1109;
  double t1110;
  double t1111;
  double t1112;
  double t1103;
  double t1104;
  double t1105;
  double t1094;
  double t1095;
  double t1096;
  double t1097;
  double t1098;
  double t1099;
  double t1122;
  double t1123;
  double t1125;
  double t1126;
  double t1127;
  double t1148;
  double t1149;
  double t1150;
  double t1144;
  double t1145;
  double t1146;
  double t1106;
  double t1113;
  double t1114;
  double t1116;
  double t1117;
  double t1118;
  double t1163;
  double t1164;
  double t1165;
  double t1166;
  double t1167;
  double t1131;
  double t1154;
  double t1155;
  double t1156;
  double t1147;
  double t1151;
  double t1152;
  double t1190;
  double t1191;
  double t1192;
  double t1162;
  double t1168;
  double t1169;
  double t1170;
  double t1171;
  double t1214;
  double t1215;
  double t1209;
  double t1210;
  double t1211;
  double t1172;
  double t1175;
  double t1176;
  double t1177;
  double t1178;
  double t1179;
  double t1180;
  double t1181;
  double t1182;
  double t1183;
  double t1184;
  double t1185;
  double t1186;
  double t1187;
  double t1188;
  double t1194;
  double t1197;
  double t1199;
  double t1229;
  double t1230;
  double t1231;
  double t1225;
  double t1226;
  double t1227;
  double t1201;
  t1068 = Cos(var1[3]);
  t1053 = Cos(var1[4]);
  t1058 = Sin(var1[3]);
  t1074 = Sin(var1[4]);
  t1087 = Cos(var1[2]);
  t1050 = Sin(var1[2]);
  t1088 = t1068*t1053;
  t1089 = -1.*t1058*t1074;
  t1090 = t1088 + t1089;
  t1066 = t1053*t1058;
  t1081 = t1068*t1074;
  t1083 = t1066 + t1081;
  t1084 = -1.*t1050*t1083;
  t1091 = t1087*t1090;
  t1092 = t1084 + t1091;
  t1107 = -1.*t1053;
  t1108 = 1. + t1107;
  t1109 = -0.16*t1108;
  t1110 = -0.167371*t1053;
  t1111 = 0.022663*t1074;
  t1112 = t1109 + t1110 + t1111;
  t1103 = -0.022663*t1053;
  t1104 = -0.007370999999999989*t1074;
  t1105 = t1103 + t1104;
  t1094 = -1.*t1053*t1058;
  t1095 = -1.*t1068*t1074;
  t1096 = t1094 + t1095;
  t1097 = t1087*t1096;
  t1098 = -1.*t1050*t1090;
  t1099 = t1097 + t1098;
  t1122 = t1050*t1096;
  t1123 = t1122 + t1091;
  t1125 = t1087*t1083;
  t1126 = t1050*t1090;
  t1127 = t1125 + t1126;
  t1148 = -1.*t1058*t1105;
  t1149 = t1068*t1112;
  t1150 = t1148 + t1149;
  t1144 = t1068*t1105;
  t1145 = t1058*t1112;
  t1146 = t1144 + t1145;
  t1106 = -1.*t1053*t1105;
  t1113 = t1112*t1074;
  t1114 = t1106 + t1113;
  t1116 = t1053*t1112;
  t1117 = t1105*t1074;
  t1118 = t1116 + t1117;
  t1163 = -1.*t1068*t1053;
  t1164 = t1058*t1074;
  t1165 = t1163 + t1164;
  t1166 = t1050*t1165;
  t1167 = t1097 + t1166;
  t1131 = -1.*t1050*t1096;
  t1154 = -1.*t1150*t1096;
  t1155 = -1.*t1146*t1090;
  t1156 = t1154 + t1155;
  t1147 = t1146*t1083;
  t1151 = t1150*t1090;
  t1152 = t1147 + t1151;
  t1190 = -1.*t1068*t1105;
  t1191 = -1.*t1058*t1112;
  t1192 = t1190 + t1191;
  t1162 = 0.0033980902199999994*t1123;
  t1168 = -0.0011052077399999983*t1167;
  t1169 = t1162 + t1168;
  t1170 = 0.5*var2[4]*t1169;
  t1171 = 0.14994*t1114*t1123;
  t1214 = -0.007370999999999989*t1053;
  t1215 = t1214 + t1111;
  t1209 = 0.022663*t1053;
  t1210 = 0.007370999999999989*t1074;
  t1211 = t1209 + t1210;
  t1172 = 0.14994*t1118*t1167;
  t1175 = 0.29988*t1123*t1127;
  t1176 = 0.29988*t1123*t1167;
  t1177 = t1175 + t1176;
  t1178 = 0.5*var2[0]*t1177;
  t1179 = 0.14994*t1123*t1092;
  t1180 = 0.14994*t1099*t1127;
  t1181 = t1087*t1165;
  t1182 = t1131 + t1181;
  t1183 = 0.14994*t1123*t1182;
  t1184 = 0.14994*t1099*t1167;
  t1185 = t1179 + t1180 + t1183 + t1184;
  t1186 = 0.5*var2[1]*t1185;
  t1187 = 0.14994*t1123*t1156;
  t1188 = t1150*t1096;
  t1194 = t1146*t1090;
  t1197 = 0.14994*t1152*t1167;
  t1199 = -1.*t1146*t1096;
  t1229 = t1068*t1211;
  t1230 = -1.*t1058*t1215;
  t1231 = t1229 + t1230;
  t1225 = t1058*t1211;
  t1226 = t1068*t1215;
  t1227 = t1225 + t1226;
  t1201 = -1.*t1150*t1165;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(0.5*(0.29988*t1099*t1123 + 0.29988*t1092*t1127)*var2[0] + 0.5*(0.14994*Power(t1092,2) + 0.14994*Power(t1099,2) + 0.14994*(-1.*t1083*t1087 + t1098)*t1127 + 0.14994*t1123*(-1.*t1087*t1090 + t1131))*var2[1] + 0.5*(0.14994*t1099*t1152 + 0.14994*t1092*t1156)*var2[2] + 0.5*(0.14994*t1092*t1114 + 0.14994*t1099*t1118)*var2[3] + 0.5*(0.0033980902199999994*t1092 - 0.0011052077399999983*t1099)*var2[4]);
  p_output1[3]=var2[0]*(t1170 + t1178 + t1186 + 0.5*(t1187 + 0.14994*t1123*(t1083*t1150 + t1188 + t1090*t1192 + t1194) + t1197 + 0.14994*t1127*(-1.*t1090*t1150 - 1.*t1096*t1192 + t1199 + t1201))*var2[2] + 0.5*(t1171 + t1172)*var2[3]);
  p_output1[4]=var2[0]*(t1170 + t1178 + t1186 + 0.5*(t1187 + t1197 + 0.14994*t1123*(t1188 + t1194 + t1083*t1227 + t1090*t1231) + 0.14994*t1127*(t1199 + t1201 - 1.*t1090*t1227 - 1.*t1096*t1231))*var2[2] + 0.5*(t1171 + t1172 + 0.14994*t1127*(t1116 + t1117 + t1074*t1211 - 1.*t1053*t1215) + 0.14994*t1123*(t1053*t1105 - 1.*t1074*t1112 + t1053*t1211 + t1074*t1215))*var2[3]);
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

#include "Ce3_vec_L4_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L4_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
