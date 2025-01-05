/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:27 GMT-05:00
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
  double t1097;
  double t1066;
  double t1076;
  double t1098;
  double t1118;
  double t1061;
  double t1119;
  double t1123;
  double t1131;
  double t1173;
  double t1175;
  double t1089;
  double t1108;
  double t1110;
  double t1115;
  double t1133;
  double t1137;
  double t1187;
  double t1188;
  double t1189;
  double t1190;
  double t1176;
  double t1179;
  double t1181;
  double t1183;
  double t1141;
  double t1145;
  double t1148;
  double t1150;
  double t1151;
  double t1156;
  double t1206;
  double t1207;
  double t1209;
  double t1210;
  double t1211;
  double t1229;
  double t1230;
  double t1237;
  double t1238;
  double t1239;
  double t1240;
  double t1241;
  double t1231;
  double t1232;
  double t1233;
  double t1234;
  double t1235;
  double t1196;
  double t1198;
  double t1199;
  double t1200;
  double t1201;
  double t1169;
  double t1184;
  double t1185;
  double t1193;
  double t1194;
  double t1254;
  double t1255;
  double t1256;
  double t1257;
  double t1258;
  double t1259;
  double t1245;
  double t1246;
  double t1247;
  double t1236;
  double t1242;
  double t1243;
  double t1283;
  double t1284;
  double t1285;
  double t1286;
  double t1287;
  double t1280;
  double t1281;
  double t1253;
  double t1260;
  double t1261;
  double t1262;
  double t1263;
  double t1307;
  double t1308;
  double t1309;
  double t1304;
  double t1305;
  double t1264;
  double t1267;
  double t1268;
  double t1269;
  double t1270;
  double t1271;
  double t1272;
  double t1273;
  double t1274;
  double t1275;
  double t1276;
  double t1277;
  double t1278;
  double t1292;
  double t1279;
  double t1289;
  double t1293;
  double t1295;
  double t1327;
  double t1328;
  double t1329;
  double t1323;
  double t1324;
  double t1325;
  double t1297;
  t1097 = Cos(var1[3]);
  t1066 = Cos(var1[4]);
  t1076 = Sin(var1[3]);
  t1098 = Sin(var1[4]);
  t1118 = Cos(var1[2]);
  t1061 = Sin(var1[2]);
  t1119 = t1097*t1066;
  t1123 = -1.*t1076*t1098;
  t1131 = t1119 + t1123;
  t1173 = -1.*t1066;
  t1175 = 1. + t1173;
  t1089 = -1.*t1066*t1076;
  t1108 = -1.*t1097*t1098;
  t1110 = t1089 + t1108;
  t1115 = -1.*t1061*t1110;
  t1133 = t1118*t1131;
  t1137 = t1115 + t1133;
  t1187 = -0.2375*t1175;
  t1188 = -0.314514*t1066;
  t1189 = 0.0012709999999999978*t1098;
  t1190 = t1187 + t1188 + t1189;
  t1176 = -0.0265*t1175;
  t1179 = -0.025229*t1066;
  t1181 = 0.07701400000000003*t1098;
  t1183 = t1176 + t1179 + t1181;
  t1141 = t1066*t1076;
  t1145 = t1097*t1098;
  t1148 = t1141 + t1145;
  t1150 = t1118*t1148;
  t1151 = -1.*t1061*t1131;
  t1156 = t1150 + t1151;
  t1206 = t1061*t1148;
  t1207 = t1206 + t1133;
  t1209 = t1118*t1110;
  t1210 = t1061*t1131;
  t1211 = t1209 + t1210;
  t1229 = -1.*t1097;
  t1230 = 1. + t1229;
  t1237 = -0.0265*t1230;
  t1238 = -0.0695*t1076;
  t1239 = -1.*t1076*t1190;
  t1240 = t1097*t1183;
  t1241 = t1237 + t1238 + t1239 + t1240;
  t1231 = -0.0695*t1230;
  t1232 = 0.0265*t1076;
  t1233 = t1097*t1190;
  t1234 = t1076*t1183;
  t1235 = t1231 + t1232 + t1233 + t1234;
  t1196 = -0.0695*t1066;
  t1198 = -1.*t1066*t1190;
  t1199 = 0.0265*t1098;
  t1200 = t1183*t1098;
  t1201 = t1196 + t1198 + t1199 + t1200;
  t1169 = 0.0265*t1066;
  t1184 = t1066*t1183;
  t1185 = 0.0695*t1098;
  t1193 = t1190*t1098;
  t1194 = t1169 + t1184 + t1185 + t1193;
  t1254 = t1061*t1110;
  t1255 = -1.*t1097*t1066;
  t1256 = t1076*t1098;
  t1257 = t1255 + t1256;
  t1258 = t1118*t1257;
  t1259 = t1254 + t1258;
  t1245 = t1241*t1110;
  t1246 = t1235*t1131;
  t1247 = t1245 + t1246;
  t1236 = -1.*t1235*t1148;
  t1242 = -1.*t1241*t1131;
  t1243 = t1236 + t1242;
  t1283 = -0.0695*t1097;
  t1284 = -0.0265*t1076;
  t1285 = -1.*t1097*t1190;
  t1286 = -1.*t1076*t1183;
  t1287 = t1283 + t1284 + t1285 + t1286;
  t1280 = 0.0265*t1097;
  t1281 = t1280 + t1238 + t1239 + t1240;
  t1253 = 0.015375074960000006*t1211;
  t1260 = 0.0002537424399999996*t1259;
  t1261 = t1253 + t1260;
  t1262 = 0.5*var2[4]*t1261;
  t1263 = 0.19964*t1201*t1211;
  t1307 = 0.07701400000000003*t1066;
  t1308 = -0.0012709999999999978*t1098;
  t1309 = t1307 + t1308;
  t1304 = 0.0012709999999999978*t1066;
  t1305 = t1304 + t1181;
  t1264 = 0.19964*t1194*t1259;
  t1267 = 0.39928*t1207*t1211;
  t1268 = 0.39928*t1211*t1259;
  t1269 = t1267 + t1268;
  t1270 = 0.5*var2[0]*t1269;
  t1271 = 0.19964*t1137*t1207;
  t1272 = 0.19964*t1156*t1211;
  t1273 = 0.19964*t1137*t1259;
  t1274 = -1.*t1061*t1257;
  t1275 = t1209 + t1274;
  t1276 = 0.19964*t1211*t1275;
  t1277 = t1271 + t1272 + t1273 + t1276;
  t1278 = 0.5*var2[1]*t1277;
  t1292 = 0.19964*t1211*t1247;
  t1279 = -1.*t1241*t1110;
  t1289 = -1.*t1235*t1131;
  t1293 = 0.19964*t1243*t1259;
  t1295 = t1235*t1110;
  t1327 = t1097*t1309;
  t1328 = -1.*t1076*t1305;
  t1329 = t1327 + t1328;
  t1323 = t1076*t1309;
  t1324 = t1097*t1305;
  t1325 = t1323 + t1324;
  t1297 = t1241*t1257;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(0.5*(0.39928*t1156*t1207 + 0.39928*t1137*t1211)*var2[0] + 0.5*(0.19964*Power(t1137,2) + 0.19964*Power(t1156,2) + 0.19964*(-1.*t1118*t1131 - 1.*t1061*t1148)*t1207 + 0.19964*(-1.*t1110*t1118 + t1151)*t1211)*var2[1] + 0.5*(0.19964*t1137*t1243 + 0.19964*t1156*t1247)*var2[2] + 0.5*(0.19964*t1137*t1194 + 0.19964*t1156*t1201)*var2[3] + 0.5*(0.0002537424399999996*t1137 + 0.015375074960000006*t1156)*var2[4]);
  p_output1[3]=var2[0]*(t1262 + t1270 + t1278 + 0.5*(0.19964*t1211*(t1279 - 1.*t1148*t1281 - 1.*t1131*t1287 + t1289) + t1292 + t1293 + 0.19964*t1207*(t1131*t1281 + t1110*t1287 + t1295 + t1297))*var2[2] + 0.5*(t1263 + t1264)*var2[3]);
  p_output1[4]=var2[0]*(t1262 + t1270 + t1278 + 0.5*(t1292 + t1293 + 0.19964*t1207*(t1295 + t1297 + t1131*t1325 + t1110*t1329) + 0.19964*t1211*(t1279 + t1289 - 1.*t1148*t1325 - 1.*t1131*t1329))*var2[2] + 0.5*(t1263 + t1264 + 0.19964*t1211*(0.0695*t1066 - 0.0265*t1098 - 1.*t1098*t1183 + t1066*t1190 + t1098*t1305 + t1066*t1309) + 0.19964*t1207*(t1169 + t1184 + t1185 + t1193 - 1.*t1066*t1305 + t1098*t1309))*var2[3]);
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

#include "Ce3_vec_L3_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L3_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
