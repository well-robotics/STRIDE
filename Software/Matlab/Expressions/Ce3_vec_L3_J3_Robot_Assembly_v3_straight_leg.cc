/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:32 GMT-05:00
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
  double t1188;
  double t1150;
  double t1167;
  double t1200;
  double t1221;
  double t1226;
  double t1216;
  double t1217;
  double t1227;
  double t1228;
  double t1231;
  double t1233;
  double t1237;
  double t1238;
  double t1240;
  double t1243;
  double t1208;
  double t1209;
  double t1210;
  double t1205;
  double t1246;
  double t1247;
  double t1248;
  double t1139;
  double t1250;
  double t1251;
  double t1252;
  double t1253;
  double t1254;
  double t1187;
  double t1201;
  double t1202;
  double t1219;
  double t1220;
  double t1234;
  double t1244;
  double t1245;
  double t1249;
  double t1256;
  double t1258;
  double t1213;
  double t1271;
  double t1272;
  double t1276;
  double t1315;
  double t1316;
  double t1319;
  double t1320;
  double t1321;
  double t1307;
  double t1313;
  double t1306;
  double t1314;
  double t1322;
  double t1323;
  double t1324;
  double t1332;
  double t1334;
  double t1335;
  double t1336;
  double t1340;
  double t1342;
  double t1343;
  double t1344;
  double t1363;
  double t1364;
  double t1366;
  double t1292;
  double t1291;
  double t1293;
  double t1295;
  double t1296;
  double t1395;
  double t1396;
  double t1397;
  double t1399;
  double t1400;
  double t1404;
  double t1405;
  double t1406;
  double t1398;
  double t1401;
  double t1402;
  double t1348;
  double t1349;
  double t1350;
  double t1351;
  double t1352;
  double t1403;
  double t1407;
  double t1408;
  double t1354;
  double t1355;
  double t1356;
  double t1357;
  double t1358;
  double t1410;
  double t1411;
  double t1412;
  double t1368;
  double t1369;
  double t1370;
  double t1371;
  double t1372;
  double t1373;
  double t1374;
  double t1382;
  double t1383;
  double t1384;
  double t1385;
  t1188 = Cos(var1[3]);
  t1150 = Cos(var1[4]);
  t1167 = Sin(var1[3]);
  t1200 = Sin(var1[4]);
  t1221 = -1.*t1150;
  t1226 = 1. + t1221;
  t1216 = -1.*t1188;
  t1217 = 1. + t1216;
  t1227 = -0.2375*t1226;
  t1228 = -0.314514*t1150;
  t1231 = 0.0012709999999999978*t1200;
  t1233 = t1227 + t1228 + t1231;
  t1237 = -0.0265*t1226;
  t1238 = -0.025229*t1150;
  t1240 = 0.07701400000000003*t1200;
  t1243 = t1237 + t1238 + t1240;
  t1208 = t1188*t1150;
  t1209 = -1.*t1167*t1200;
  t1210 = t1208 + t1209;
  t1205 = Sin(var1[2]);
  t1246 = t1150*t1167;
  t1247 = t1188*t1200;
  t1248 = t1246 + t1247;
  t1139 = Cos(var1[2]);
  t1250 = -0.0265*t1217;
  t1251 = -0.0695*t1167;
  t1252 = -1.*t1167*t1233;
  t1253 = t1188*t1243;
  t1254 = t1250 + t1251 + t1252 + t1253;
  t1187 = -1.*t1150*t1167;
  t1201 = -1.*t1188*t1200;
  t1202 = t1187 + t1201;
  t1219 = -0.0695*t1217;
  t1220 = 0.0265*t1167;
  t1234 = t1188*t1233;
  t1244 = t1167*t1243;
  t1245 = t1219 + t1220 + t1234 + t1244;
  t1249 = -1.*t1245*t1248;
  t1256 = -1.*t1254*t1210;
  t1258 = t1249 + t1256;
  t1213 = -1.*t1205*t1210;
  t1271 = t1254*t1202;
  t1272 = t1245*t1210;
  t1276 = t1271 + t1272;
  t1315 = -0.0695*t1188;
  t1316 = -0.0265*t1167;
  t1319 = -1.*t1188*t1233;
  t1320 = -1.*t1167*t1243;
  t1321 = t1315 + t1316 + t1319 + t1320;
  t1307 = 0.0265*t1188;
  t1313 = t1307 + t1251 + t1252 + t1253;
  t1306 = -1.*t1254*t1202;
  t1314 = -1.*t1313*t1248;
  t1322 = -1.*t1321*t1210;
  t1323 = -1.*t1245*t1210;
  t1324 = t1306 + t1314 + t1322 + t1323;
  t1332 = t1321*t1202;
  t1334 = t1245*t1202;
  t1335 = t1313*t1210;
  t1336 = -1.*t1188*t1150;
  t1340 = t1167*t1200;
  t1342 = t1336 + t1340;
  t1343 = t1254*t1342;
  t1344 = t1332 + t1334 + t1335 + t1343;
  t1363 = t1139*t1202;
  t1364 = t1205*t1210;
  t1366 = t1363 + t1364;
  t1292 = t1139*t1210;
  t1291 = -1.*t1205*t1202;
  t1293 = t1291 + t1292;
  t1295 = t1139*t1248;
  t1296 = t1295 + t1213;
  t1395 = 0.07701400000000003*t1150;
  t1396 = -0.0012709999999999978*t1200;
  t1397 = t1395 + t1396;
  t1399 = 0.0012709999999999978*t1150;
  t1400 = t1399 + t1240;
  t1404 = t1188*t1397;
  t1405 = -1.*t1167*t1400;
  t1406 = t1404 + t1405;
  t1398 = t1167*t1397;
  t1401 = t1188*t1400;
  t1402 = t1398 + t1401;
  t1348 = 0.0265*t1150;
  t1349 = t1150*t1243;
  t1350 = 0.0695*t1200;
  t1351 = t1233*t1200;
  t1352 = t1348 + t1349 + t1350 + t1351;
  t1403 = -1.*t1402*t1248;
  t1407 = -1.*t1406*t1210;
  t1408 = t1306 + t1403 + t1323 + t1407;
  t1354 = -0.0695*t1150;
  t1355 = -1.*t1150*t1233;
  t1356 = 0.0265*t1200;
  t1357 = t1243*t1200;
  t1358 = t1354 + t1355 + t1356 + t1357;
  t1410 = t1406*t1202;
  t1411 = t1402*t1210;
  t1412 = t1334 + t1410 + t1411 + t1343;
  t1368 = 0.19964*t1366*t1276;
  t1369 = t1205*t1202;
  t1370 = t1139*t1342;
  t1371 = t1369 + t1370;
  t1372 = 0.19964*t1258*t1371;
  t1373 = t1205*t1248;
  t1374 = t1373 + t1292;
  t1382 = 0.19964*t1293*t1276;
  t1383 = -1.*t1205*t1342;
  t1384 = t1363 + t1383;
  t1385 = 0.19964*t1258*t1384;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.19964*t1258*t1293 + 0.19964*t1276*t1296)*var2[0] + 0.5*(0.19964*(-1.*t1139*t1202 + t1213)*t1258 + 0.19964*(-1.*t1139*t1210 - 1.*t1205*t1248)*t1276)*var2[1])*var2[2];
  p_output1[3]=var2[2]*(0.5*(0.19964*t1324*t1366 + t1368 + t1372 + 0.19964*t1344*t1374)*var2[0] + 0.5*(0.19964*t1293*t1324 + 0.19964*t1296*t1344 + t1382 + t1385)*var2[1] + 0.5*(0.39928*t1258*t1324 + 0.39928*t1276*t1344)*var2[2] + 0.5*(0.19964*t1324*t1352 + 0.19964*t1344*t1358)*var2[3] + 0.5*(0.0002537424399999996*t1324 + 0.015375074960000006*t1344)*var2[4]);
  p_output1[4]=var2[2]*(0.5*(t1368 + t1372 + 0.19964*t1366*t1408 + 0.19964*t1374*t1412)*var2[0] + 0.5*(t1382 + t1385 + 0.19964*t1293*t1408 + 0.19964*t1296*t1412)*var2[1] + 0.5*(0.39928*t1258*t1408 + 0.39928*t1276*t1412)*var2[2] + 0.5*(0.19964*t1276*(t1348 + t1349 + t1350 + t1351 + t1200*t1397 - 1.*t1150*t1400) + 0.19964*t1258*(0.0695*t1150 - 0.0265*t1200 + t1150*t1233 - 1.*t1200*t1243 + t1150*t1397 + t1200*t1400) + 0.19964*t1352*t1408 + 0.19964*t1358*t1412)*var2[3] + 0.5*(0.0002537424399999996*t1408 + 0.015375074960000006*t1412)*var2[4]);
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

#include "Ce3_vec_L3_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L3_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
