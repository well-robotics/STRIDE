/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:29 GMT-05:00
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
  double t1141;
  double t1133;
  double t1137;
  double t1145;
  double t1160;
  double t1119;
  double t1167;
  double t1176;
  double t1179;
  double t1206;
  double t1207;
  double t1139;
  double t1150;
  double t1156;
  double t1157;
  double t1181;
  double t1187;
  double t1214;
  double t1216;
  double t1217;
  double t1218;
  double t1208;
  double t1209;
  double t1210;
  double t1211;
  double t1189;
  double t1194;
  double t1195;
  double t1196;
  double t1199;
  double t1200;
  double t1232;
  double t1233;
  double t1234;
  double t1238;
  double t1240;
  double t1265;
  double t1266;
  double t1276;
  double t1280;
  double t1282;
  double t1283;
  double t1284;
  double t1267;
  double t1268;
  double t1271;
  double t1272;
  double t1273;
  double t1253;
  double t1205;
  double t1212;
  double t1213;
  double t1219;
  double t1220;
  double t1222;
  double t1223;
  double t1224;
  double t1225;
  double t1226;
  double t1301;
  double t1302;
  double t1303;
  double t1304;
  double t1306;
  double t1248;
  double t1249;
  double t1254;
  double t1256;
  double t1292;
  double t1293;
  double t1294;
  double t1275;
  double t1288;
  double t1290;
  double t1335;
  double t1336;
  double t1337;
  double t1338;
  double t1339;
  double t1332;
  double t1333;
  double t1300;
  double t1307;
  double t1308;
  double t1310;
  double t1311;
  double t1364;
  double t1365;
  double t1357;
  double t1358;
  double t1359;
  double t1312;
  double t1315;
  double t1316;
  double t1317;
  double t1318;
  double t1319;
  double t1320;
  double t1321;
  double t1322;
  double t1323;
  double t1324;
  double t1326;
  double t1327;
  double t1330;
  double t1344;
  double t1331;
  double t1341;
  double t1345;
  double t1347;
  double t1379;
  double t1380;
  double t1381;
  double t1375;
  double t1376;
  double t1377;
  double t1349;
  t1141 = Cos(var1[3]);
  t1133 = Cos(var1[4]);
  t1137 = Sin(var1[3]);
  t1145 = Sin(var1[4]);
  t1160 = Cos(var1[2]);
  t1119 = Sin(var1[2]);
  t1167 = t1141*t1133;
  t1176 = -1.*t1137*t1145;
  t1179 = t1167 + t1176;
  t1206 = -1.*t1133;
  t1207 = 1. + t1206;
  t1139 = t1133*t1137;
  t1150 = t1141*t1145;
  t1156 = t1139 + t1150;
  t1157 = -1.*t1119*t1156;
  t1181 = -1.*t1160*t1179;
  t1187 = t1157 + t1181;
  t1214 = -0.0265*t1207;
  t1216 = -0.025229*t1133;
  t1217 = 0.07701400000000003*t1145;
  t1218 = t1214 + t1216 + t1217;
  t1208 = -0.2375*t1207;
  t1209 = -0.314514*t1133;
  t1210 = 0.0012709999999999978*t1145;
  t1211 = t1208 + t1209 + t1210;
  t1189 = -1.*t1133*t1137;
  t1194 = -1.*t1141*t1145;
  t1195 = t1189 + t1194;
  t1196 = -1.*t1160*t1195;
  t1199 = -1.*t1119*t1179;
  t1200 = t1196 + t1199;
  t1232 = -1.*t1119*t1195;
  t1233 = t1160*t1179;
  t1234 = t1232 + t1233;
  t1238 = t1160*t1156;
  t1240 = t1238 + t1199;
  t1265 = -1.*t1141;
  t1266 = 1. + t1265;
  t1276 = -0.0265*t1266;
  t1280 = -0.0695*t1137;
  t1282 = -1.*t1137*t1211;
  t1283 = t1141*t1218;
  t1284 = t1276 + t1280 + t1282 + t1283;
  t1267 = -0.0695*t1266;
  t1268 = 0.0265*t1137;
  t1271 = t1141*t1211;
  t1272 = t1137*t1218;
  t1273 = t1267 + t1268 + t1271 + t1272;
  t1253 = t1160*t1195;
  t1205 = -0.0695*t1133;
  t1212 = -1.*t1133*t1211;
  t1213 = 0.0265*t1145;
  t1219 = t1218*t1145;
  t1220 = t1205 + t1212 + t1213 + t1219;
  t1222 = 0.0265*t1133;
  t1223 = t1133*t1218;
  t1224 = 0.0695*t1145;
  t1225 = t1211*t1145;
  t1226 = t1222 + t1223 + t1224 + t1225;
  t1301 = -1.*t1141*t1133;
  t1302 = t1137*t1145;
  t1303 = t1301 + t1302;
  t1304 = -1.*t1119*t1303;
  t1306 = t1253 + t1304;
  t1248 = t1119*t1156;
  t1249 = t1248 + t1233;
  t1254 = t1119*t1179;
  t1256 = t1253 + t1254;
  t1292 = t1284*t1195;
  t1293 = t1273*t1179;
  t1294 = t1292 + t1293;
  t1275 = -1.*t1273*t1156;
  t1288 = -1.*t1284*t1179;
  t1290 = t1275 + t1288;
  t1335 = -0.0695*t1141;
  t1336 = -0.0265*t1137;
  t1337 = -1.*t1141*t1211;
  t1338 = -1.*t1137*t1218;
  t1339 = t1335 + t1336 + t1337 + t1338;
  t1332 = 0.0265*t1141;
  t1333 = t1332 + t1280 + t1282 + t1283;
  t1300 = 0.015375074960000006*t1234;
  t1307 = 0.0002537424399999996*t1306;
  t1308 = t1300 + t1307;
  t1310 = 0.5*var2[4]*t1308;
  t1311 = 0.19964*t1220*t1234;
  t1364 = 0.0012709999999999978*t1133;
  t1365 = t1364 + t1217;
  t1357 = 0.07701400000000003*t1133;
  t1358 = -0.0012709999999999978*t1145;
  t1359 = t1357 + t1358;
  t1312 = 0.19964*t1226*t1306;
  t1315 = 0.39928*t1234*t1240;
  t1316 = 0.39928*t1234*t1306;
  t1317 = t1315 + t1316;
  t1318 = 0.5*var2[1]*t1317;
  t1319 = 0.19964*t1234*t1249;
  t1320 = 0.19964*t1240*t1256;
  t1321 = t1119*t1195;
  t1322 = t1160*t1303;
  t1323 = t1321 + t1322;
  t1324 = 0.19964*t1234*t1323;
  t1326 = 0.19964*t1256*t1306;
  t1327 = t1319 + t1320 + t1324 + t1326;
  t1330 = 0.5*var2[0]*t1327;
  t1344 = 0.19964*t1234*t1294;
  t1331 = -1.*t1284*t1195;
  t1341 = -1.*t1273*t1179;
  t1345 = 0.19964*t1290*t1306;
  t1347 = t1273*t1195;
  t1379 = t1141*t1359;
  t1380 = -1.*t1137*t1365;
  t1381 = t1379 + t1380;
  t1375 = t1137*t1359;
  t1376 = t1141*t1365;
  t1377 = t1375 + t1376;
  t1349 = t1284*t1303;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(0.5*(0.19964*Power(t1234,2) + 0.19964*Power(t1240,2) + 0.19964*t1187*t1249 + 0.19964*t1200*t1256)*var2[0] + 0.5*(0.39928*t1200*t1234 + 0.39928*t1187*t1240)*var2[1] + 0.5*(0.19964*t1200*t1290 + 0.19964*t1187*t1294)*var2[2] + 0.5*(0.19964*t1187*t1220 + 0.19964*t1200*t1226)*var2[3] + 0.5*(0.015375074960000006*t1187 + 0.0002537424399999996*t1200)*var2[4]);
  p_output1[3]=var2[1]*(t1310 + t1318 + t1330 + 0.5*(0.19964*t1234*(t1331 - 1.*t1156*t1333 - 1.*t1179*t1339 + t1341) + t1344 + t1345 + 0.19964*t1240*(t1179*t1333 + t1195*t1339 + t1347 + t1349))*var2[2] + 0.5*(t1311 + t1312)*var2[3]);
  p_output1[4]=var2[1]*(t1310 + t1318 + t1330 + 0.5*(t1344 + t1345 + 0.19964*t1234*(t1331 + t1341 - 1.*t1156*t1377 - 1.*t1179*t1381) + 0.19964*t1240*(t1347 + t1349 + t1179*t1377 + t1195*t1381))*var2[2] + 0.5*(t1311 + t1312 + 0.19964*t1240*(t1222 + t1223 + t1224 + t1225 + t1145*t1359 - 1.*t1133*t1365) + 0.19964*t1234*(0.0695*t1133 - 0.0265*t1145 + t1133*t1211 - 1.*t1145*t1218 + t1133*t1359 + t1145*t1365))*var2[3]);
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

#include "Ce3_vec_L3_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L3_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
