/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:51 GMT-05:00
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
  double t1092;
  double t1081;
  double t1088;
  double t1093;
  double t1101;
  double t1066;
  double t1102;
  double t1103;
  double t1104;
  double t1091;
  double t1097;
  double t1099;
  double t1100;
  double t1109;
  double t1110;
  double t1132;
  double t1133;
  double t1134;
  double t1124;
  double t1125;
  double t1126;
  double t1127;
  double t1128;
  double t1129;
  double t1113;
  double t1114;
  double t1115;
  double t1118;
  double t1119;
  double t1120;
  double t1151;
  double t1149;
  double t1152;
  double t1144;
  double t1145;
  double t1183;
  double t1184;
  double t1187;
  double t1179;
  double t1180;
  double t1181;
  double t1138;
  double t1139;
  double t1140;
  double t1130;
  double t1135;
  double t1136;
  double t1205;
  double t1206;
  double t1207;
  double t1208;
  double t1209;
  double t1158;
  double t1159;
  double t1167;
  double t1168;
  double t1173;
  double t1194;
  double t1195;
  double t1196;
  double t1182;
  double t1188;
  double t1189;
  double t1237;
  double t1238;
  double t1239;
  double t1204;
  double t1210;
  double t1212;
  double t1213;
  double t1214;
  double t1258;
  double t1259;
  double t1260;
  double t1255;
  double t1256;
  double t1216;
  double t1219;
  double t1220;
  double t1221;
  double t1222;
  double t1223;
  double t1224;
  double t1225;
  double t1226;
  double t1228;
  double t1229;
  double t1232;
  double t1233;
  double t1234;
  double t1235;
  double t1241;
  double t1244;
  double t1246;
  double t1276;
  double t1277;
  double t1278;
  double t1272;
  double t1273;
  double t1274;
  double t1248;
  t1092 = Cos(var1[3]);
  t1081 = Cos(var1[4]);
  t1088 = Sin(var1[3]);
  t1093 = Sin(var1[4]);
  t1101 = Cos(var1[2]);
  t1066 = Sin(var1[2]);
  t1102 = t1092*t1081;
  t1103 = -1.*t1088*t1093;
  t1104 = t1102 + t1103;
  t1091 = -1.*t1081*t1088;
  t1097 = -1.*t1092*t1093;
  t1099 = t1091 + t1097;
  t1100 = -1.*t1066*t1099;
  t1109 = -1.*t1101*t1104;
  t1110 = t1100 + t1109;
  t1132 = -0.022663*t1081;
  t1133 = -0.007370999999999989*t1093;
  t1134 = t1132 + t1133;
  t1124 = -1.*t1081;
  t1125 = 1. + t1124;
  t1126 = -0.16*t1125;
  t1127 = -0.167371*t1081;
  t1128 = 0.022663*t1093;
  t1129 = t1126 + t1127 + t1128;
  t1113 = t1081*t1088;
  t1114 = t1092*t1093;
  t1115 = t1113 + t1114;
  t1118 = -1.*t1101*t1115;
  t1119 = -1.*t1066*t1104;
  t1120 = t1118 + t1119;
  t1151 = t1101*t1104;
  t1149 = -1.*t1066*t1115;
  t1152 = t1149 + t1151;
  t1144 = t1101*t1099;
  t1145 = t1144 + t1119;
  t1183 = -1.*t1088*t1134;
  t1184 = t1092*t1129;
  t1187 = t1183 + t1184;
  t1179 = t1092*t1134;
  t1180 = t1088*t1129;
  t1181 = t1179 + t1180;
  t1138 = -1.*t1081*t1134;
  t1139 = t1129*t1093;
  t1140 = t1138 + t1139;
  t1130 = t1081*t1129;
  t1135 = t1134*t1093;
  t1136 = t1130 + t1135;
  t1205 = -1.*t1092*t1081;
  t1206 = t1088*t1093;
  t1207 = t1205 + t1206;
  t1208 = t1101*t1207;
  t1209 = t1100 + t1208;
  t1158 = t1066*t1099;
  t1159 = t1158 + t1151;
  t1167 = t1101*t1115;
  t1168 = t1066*t1104;
  t1173 = t1167 + t1168;
  t1194 = -1.*t1187*t1099;
  t1195 = -1.*t1181*t1104;
  t1196 = t1194 + t1195;
  t1182 = t1181*t1115;
  t1188 = t1187*t1104;
  t1189 = t1182 + t1188;
  t1237 = -1.*t1092*t1134;
  t1238 = -1.*t1088*t1129;
  t1239 = t1237 + t1238;
  t1204 = 0.0033980902199999994*t1145;
  t1210 = -0.0011052077399999983*t1209;
  t1212 = t1204 + t1210;
  t1213 = 0.5*var2[4]*t1212;
  t1214 = 0.14994*t1140*t1145;
  t1258 = 0.022663*t1081;
  t1259 = 0.007370999999999989*t1093;
  t1260 = t1258 + t1259;
  t1255 = -0.007370999999999989*t1081;
  t1256 = t1255 + t1128;
  t1216 = 0.14994*t1136*t1209;
  t1219 = 0.29988*t1152*t1145;
  t1220 = 0.29988*t1145*t1209;
  t1221 = t1219 + t1220;
  t1222 = 0.5*var2[1]*t1221;
  t1223 = 0.14994*t1159*t1152;
  t1224 = 0.14994*t1145*t1173;
  t1225 = 0.14994*t1159*t1209;
  t1226 = t1066*t1207;
  t1228 = t1144 + t1226;
  t1229 = 0.14994*t1145*t1228;
  t1232 = t1223 + t1224 + t1225 + t1229;
  t1233 = 0.5*var2[0]*t1232;
  t1234 = 0.14994*t1145*t1196;
  t1235 = t1187*t1099;
  t1241 = t1181*t1104;
  t1244 = 0.14994*t1189*t1209;
  t1246 = -1.*t1181*t1099;
  t1276 = t1092*t1260;
  t1277 = -1.*t1088*t1256;
  t1278 = t1276 + t1277;
  t1272 = t1088*t1260;
  t1273 = t1092*t1256;
  t1274 = t1272 + t1273;
  t1248 = -1.*t1187*t1207;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(0.5*(0.14994*Power(t1145,2) + 0.14994*Power(t1152,2) + 0.14994*t1110*t1159 + 0.14994*t1120*t1173)*var2[0] + 0.5*(0.29988*t1110*t1145 + 0.29988*t1120*t1152)*var2[1] + 0.5*(0.14994*t1110*t1189 + 0.14994*t1120*t1196)*var2[2] + 0.5*(0.14994*t1110*t1136 + 0.14994*t1120*t1140)*var2[3] + 0.5*(-0.0011052077399999983*t1110 + 0.0033980902199999994*t1120)*var2[4]);
  p_output1[3]=var2[1]*(t1213 + t1222 + t1233 + 0.5*(t1234 + 0.14994*t1145*(t1115*t1187 + t1235 + t1104*t1239 + t1241) + t1244 + 0.14994*t1152*(-1.*t1104*t1187 - 1.*t1099*t1239 + t1246 + t1248))*var2[2] + 0.5*(t1214 + t1216)*var2[3]);
  p_output1[4]=var2[1]*(t1213 + t1222 + t1233 + 0.5*(t1234 + t1244 + 0.14994*t1152*(t1246 + t1248 - 1.*t1104*t1274 - 1.*t1099*t1278) + 0.14994*t1145*(t1235 + t1241 + t1115*t1274 + t1104*t1278))*var2[2] + 0.5*(t1214 + t1216 + 0.14994*t1145*(-1.*t1093*t1129 + t1081*t1134 + t1093*t1256 + t1081*t1260) + 0.14994*t1152*(t1130 + t1135 - 1.*t1081*t1256 + t1093*t1260))*var2[3]);
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

#include "Ce3_vec_L4_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L4_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
