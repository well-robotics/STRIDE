/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:41:03 GMT-05:00
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
  double t1202;
  double t1161;
  double t1174;
  double t1203;
  double t1243;
  double t1152;
  double t1247;
  double t1251;
  double t1252;
  double t1196;
  double t1224;
  double t1241;
  double t1242;
  double t1253;
  double t1254;
  double t1291;
  double t1292;
  double t1294;
  double t1302;
  double t1305;
  double t1310;
  double t1279;
  double t1286;
  double t1287;
  double t1259;
  double t1262;
  double t1263;
  double t1265;
  double t1266;
  double t1267;
  double t1329;
  double t1330;
  double t1333;
  double t1335;
  double t1336;
  double t1362;
  double t1363;
  double t1364;
  double t1358;
  double t1359;
  double t1360;
  double t1288;
  double t1313;
  double t1314;
  double t1318;
  double t1319;
  double t1324;
  double t1380;
  double t1381;
  double t1382;
  double t1383;
  double t1384;
  double t1340;
  double t1371;
  double t1372;
  double t1373;
  double t1361;
  double t1366;
  double t1367;
  double t1407;
  double t1408;
  double t1409;
  double t1379;
  double t1385;
  double t1386;
  double t1387;
  double t1388;
  double t1431;
  double t1432;
  double t1426;
  double t1427;
  double t1428;
  double t1389;
  double t1392;
  double t1393;
  double t1394;
  double t1395;
  double t1396;
  double t1397;
  double t1398;
  double t1399;
  double t1400;
  double t1401;
  double t1402;
  double t1403;
  double t1404;
  double t1405;
  double t1411;
  double t1414;
  double t1416;
  double t1446;
  double t1447;
  double t1448;
  double t1442;
  double t1443;
  double t1444;
  double t1418;
  t1202 = Cos(var1[5]);
  t1161 = Cos(var1[6]);
  t1174 = Sin(var1[5]);
  t1203 = Sin(var1[6]);
  t1243 = Cos(var1[2]);
  t1152 = Sin(var1[2]);
  t1247 = t1202*t1161;
  t1251 = -1.*t1174*t1203;
  t1252 = t1247 + t1251;
  t1196 = t1161*t1174;
  t1224 = t1202*t1203;
  t1241 = t1196 + t1224;
  t1242 = -1.*t1152*t1241;
  t1253 = t1243*t1252;
  t1254 = t1242 + t1253;
  t1291 = -1.*t1161;
  t1292 = 1. + t1291;
  t1294 = -0.16*t1292;
  t1302 = -0.167368*t1161;
  t1305 = 0.022659*t1203;
  t1310 = t1294 + t1302 + t1305;
  t1279 = -0.022659*t1161;
  t1286 = -0.007367999999999986*t1203;
  t1287 = t1279 + t1286;
  t1259 = -1.*t1161*t1174;
  t1262 = -1.*t1202*t1203;
  t1263 = t1259 + t1262;
  t1265 = t1243*t1263;
  t1266 = -1.*t1152*t1252;
  t1267 = t1265 + t1266;
  t1329 = t1152*t1263;
  t1330 = t1329 + t1253;
  t1333 = t1243*t1241;
  t1335 = t1152*t1252;
  t1336 = t1333 + t1335;
  t1362 = -1.*t1174*t1287;
  t1363 = t1202*t1310;
  t1364 = t1362 + t1363;
  t1358 = t1202*t1287;
  t1359 = t1174*t1310;
  t1360 = t1358 + t1359;
  t1288 = -1.*t1161*t1287;
  t1313 = t1310*t1203;
  t1314 = t1288 + t1313;
  t1318 = t1161*t1310;
  t1319 = t1287*t1203;
  t1324 = t1318 + t1319;
  t1380 = -1.*t1202*t1161;
  t1381 = t1174*t1203;
  t1382 = t1380 + t1381;
  t1383 = t1152*t1382;
  t1384 = t1265 + t1383;
  t1340 = -1.*t1152*t1263;
  t1371 = -1.*t1364*t1263;
  t1372 = -1.*t1360*t1252;
  t1373 = t1371 + t1372;
  t1361 = t1360*t1241;
  t1366 = t1364*t1252;
  t1367 = t1361 + t1366;
  t1407 = -1.*t1202*t1287;
  t1408 = -1.*t1174*t1310;
  t1409 = t1407 + t1408;
  t1379 = 0.0033974904599999994*t1330;
  t1385 = -0.0011047579199999977*t1384;
  t1386 = t1379 + t1385;
  t1387 = 0.5*var2[6]*t1386;
  t1388 = 0.14994*t1314*t1330;
  t1431 = -0.007367999999999986*t1161;
  t1432 = t1431 + t1305;
  t1426 = 0.022659*t1161;
  t1427 = 0.007367999999999986*t1203;
  t1428 = t1426 + t1427;
  t1389 = 0.14994*t1324*t1384;
  t1392 = 0.29988*t1330*t1336;
  t1393 = 0.29988*t1330*t1384;
  t1394 = t1392 + t1393;
  t1395 = 0.5*var2[0]*t1394;
  t1396 = 0.14994*t1330*t1254;
  t1397 = 0.14994*t1267*t1336;
  t1398 = t1243*t1382;
  t1399 = t1340 + t1398;
  t1400 = 0.14994*t1330*t1399;
  t1401 = 0.14994*t1267*t1384;
  t1402 = t1396 + t1397 + t1400 + t1401;
  t1403 = 0.5*var2[1]*t1402;
  t1404 = 0.14994*t1330*t1373;
  t1405 = t1364*t1263;
  t1411 = t1360*t1252;
  t1414 = 0.14994*t1367*t1384;
  t1416 = -1.*t1360*t1263;
  t1446 = t1202*t1428;
  t1447 = -1.*t1174*t1432;
  t1448 = t1446 + t1447;
  t1442 = t1174*t1428;
  t1443 = t1202*t1432;
  t1444 = t1442 + t1443;
  t1418 = -1.*t1364*t1382;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(0.5*(0.29988*t1267*t1330 + 0.29988*t1254*t1336)*var2[0] + 0.5*(0.14994*Power(t1254,2) + 0.14994*Power(t1267,2) + 0.14994*(-1.*t1241*t1243 + t1266)*t1336 + 0.14994*t1330*(-1.*t1243*t1252 + t1340))*var2[1] + 0.5*(0.14994*t1267*t1367 + 0.14994*t1254*t1373)*var2[2] + 0.5*(0.14994*t1254*t1314 + 0.14994*t1267*t1324)*var2[5] + 0.5*(0.0033974904599999994*t1254 - 0.0011047579199999977*t1267)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(t1387 + t1395 + t1403 + 0.5*(t1404 + 0.14994*t1330*(t1241*t1364 + t1405 + t1252*t1409 + t1411) + t1414 + 0.14994*t1336*(-1.*t1252*t1364 - 1.*t1263*t1409 + t1416 + t1418))*var2[2] + 0.5*(t1388 + t1389)*var2[5]);
  p_output1[6]=var2[0]*(t1387 + t1395 + t1403 + 0.5*(t1404 + t1414 + 0.14994*t1330*(t1405 + t1411 + t1241*t1444 + t1252*t1448) + 0.14994*t1336*(t1416 + t1418 - 1.*t1252*t1444 - 1.*t1263*t1448))*var2[2] + 0.5*(t1388 + t1389 + 0.14994*t1336*(t1318 + t1319 + t1203*t1428 - 1.*t1161*t1432) + 0.14994*t1330*(t1161*t1287 - 1.*t1203*t1310 + t1161*t1428 + t1203*t1432))*var2[5]);
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

#include "Ce3_vec_L5_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L5_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
