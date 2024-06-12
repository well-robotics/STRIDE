/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:53 GMT-05:00
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
  double t1114;
  double t1110;
  double t1111;
  double t1120;
  double t1137;
  double t1139;
  double t1140;
  double t1142;
  double t1143;
  double t1144;
  double t1145;
  double t1147;
  double t1151;
  double t1127;
  double t1128;
  double t1132;
  double t1126;
  double t1156;
  double t1157;
  double t1158;
  double t1102;
  double t1160;
  double t1161;
  double t1162;
  double t1113;
  double t1121;
  double t1122;
  double t1141;
  double t1152;
  double t1153;
  double t1159;
  double t1164;
  double t1166;
  double t1133;
  double t1175;
  double t1176;
  double t1179;
  double t1210;
  double t1217;
  double t1218;
  double t1208;
  double t1209;
  double t1219;
  double t1220;
  double t1223;
  double t1225;
  double t1226;
  double t1228;
  double t1229;
  double t1234;
  double t1235;
  double t1236;
  double t1240;
  double t1123;
  double t1136;
  double t1188;
  double t1168;
  double t1173;
  double t1174;
  double t1266;
  double t1267;
  double t1288;
  double t1289;
  double t1290;
  double t1292;
  double t1293;
  double t1297;
  double t1298;
  double t1299;
  double t1291;
  double t1294;
  double t1295;
  double t1244;
  double t1245;
  double t1247;
  double t1296;
  double t1300;
  double t1301;
  double t1250;
  double t1251;
  double t1252;
  double t1303;
  double t1304;
  double t1305;
  double t1257;
  double t1259;
  double t1261;
  double t1262;
  double t1268;
  double t1270;
  double t1271;
  double t1272;
  double t1273;
  double t1275;
  double t1276;
  t1114 = Cos(var1[3]);
  t1110 = Cos(var1[4]);
  t1111 = Sin(var1[3]);
  t1120 = Sin(var1[4]);
  t1137 = -0.022663*t1110;
  t1139 = -0.007370999999999989*t1120;
  t1140 = t1137 + t1139;
  t1142 = -1.*t1110;
  t1143 = 1. + t1142;
  t1144 = -0.16*t1143;
  t1145 = -0.167371*t1110;
  t1147 = 0.022663*t1120;
  t1151 = t1144 + t1145 + t1147;
  t1127 = t1114*t1110;
  t1128 = -1.*t1111*t1120;
  t1132 = t1127 + t1128;
  t1126 = Sin(var1[2]);
  t1156 = t1110*t1111;
  t1157 = t1114*t1120;
  t1158 = t1156 + t1157;
  t1102 = Cos(var1[2]);
  t1160 = -1.*t1111*t1140;
  t1161 = t1114*t1151;
  t1162 = t1160 + t1161;
  t1113 = -1.*t1110*t1111;
  t1121 = -1.*t1114*t1120;
  t1122 = t1113 + t1121;
  t1141 = t1114*t1140;
  t1152 = t1111*t1151;
  t1153 = t1141 + t1152;
  t1159 = t1153*t1158;
  t1164 = t1162*t1132;
  t1166 = t1159 + t1164;
  t1133 = -1.*t1126*t1132;
  t1175 = -1.*t1162*t1122;
  t1176 = -1.*t1153*t1132;
  t1179 = t1175 + t1176;
  t1210 = -1.*t1114*t1140;
  t1217 = -1.*t1111*t1151;
  t1218 = t1210 + t1217;
  t1208 = t1162*t1122;
  t1209 = t1162*t1158;
  t1219 = t1218*t1132;
  t1220 = t1153*t1132;
  t1223 = t1208 + t1209 + t1219 + t1220;
  t1225 = -1.*t1218*t1122;
  t1226 = -1.*t1153*t1122;
  t1228 = -1.*t1162*t1132;
  t1229 = -1.*t1114*t1110;
  t1234 = t1111*t1120;
  t1235 = t1229 + t1234;
  t1236 = -1.*t1162*t1235;
  t1240 = t1225 + t1226 + t1228 + t1236;
  t1123 = t1102*t1122;
  t1136 = t1123 + t1133;
  t1188 = -1.*t1126*t1122;
  t1168 = -1.*t1126*t1158;
  t1173 = t1102*t1132;
  t1174 = t1168 + t1173;
  t1266 = t1126*t1122;
  t1267 = t1266 + t1173;
  t1288 = 0.022663*t1110;
  t1289 = 0.007370999999999989*t1120;
  t1290 = t1288 + t1289;
  t1292 = -0.007370999999999989*t1110;
  t1293 = t1292 + t1147;
  t1297 = t1114*t1290;
  t1298 = -1.*t1111*t1293;
  t1299 = t1297 + t1298;
  t1291 = t1111*t1290;
  t1294 = t1114*t1293;
  t1295 = t1291 + t1294;
  t1244 = t1110*t1151;
  t1245 = t1140*t1120;
  t1247 = t1244 + t1245;
  t1296 = t1295*t1158;
  t1300 = t1299*t1132;
  t1301 = t1208 + t1296 + t1220 + t1300;
  t1250 = -1.*t1110*t1140;
  t1251 = t1151*t1120;
  t1252 = t1250 + t1251;
  t1303 = -1.*t1299*t1122;
  t1304 = -1.*t1295*t1132;
  t1305 = t1226 + t1303 + t1304 + t1236;
  t1257 = 0.14994*t1136*t1179;
  t1259 = t1102*t1235;
  t1261 = t1188 + t1259;
  t1262 = 0.14994*t1166*t1261;
  t1268 = 0.14994*t1267*t1179;
  t1270 = t1126*t1235;
  t1271 = t1123 + t1270;
  t1272 = 0.14994*t1166*t1271;
  t1273 = t1102*t1158;
  t1275 = t1126*t1132;
  t1276 = t1273 + t1275;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.14994*t1136*t1166 + 0.14994*t1174*t1179)*var2[0] + 0.5*(0.14994*(t1133 - 1.*t1102*t1158)*t1179 + 0.14994*t1166*(-1.*t1102*t1132 + t1188))*var2[1])*var2[2];
  p_output1[3]=var2[2]*(0.5*(0.14994*t1223*t1267 + t1268 + t1272 + 0.14994*t1240*t1276)*var2[0] + 0.5*(0.14994*t1136*t1223 + 0.14994*t1174*t1240 + t1257 + t1262)*var2[1] + 0.5*(0.29988*t1166*t1223 + 0.29988*t1179*t1240)*var2[2] + 0.5*(0.14994*t1223*t1247 + 0.14994*t1240*t1252)*var2[3] + 0.5*(-0.0011052077399999983*t1223 + 0.0033980902199999994*t1240)*var2[4]);
  p_output1[4]=var2[2]*(0.5*(t1268 + t1272 + 0.14994*t1267*t1301 + 0.14994*t1276*t1305)*var2[0] + 0.5*(t1257 + t1262 + 0.14994*t1136*t1301 + 0.14994*t1174*t1305)*var2[1] + 0.5*(0.29988*t1166*t1301 + 0.29988*t1179*t1305)*var2[2] + 0.5*(0.14994*t1179*(t1244 + t1245 + t1120*t1290 - 1.*t1110*t1293) + 0.14994*t1166*(t1110*t1140 - 1.*t1120*t1151 + t1110*t1290 + t1120*t1293) + 0.14994*t1247*t1301 + 0.14994*t1252*t1305)*var2[3] + 0.5*(-0.0011052077399999983*t1301 + 0.0033980902199999994*t1305)*var2[4]);
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

#include "Ce3_vec_L4_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L4_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
