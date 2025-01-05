/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:34 GMT-05:00
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
  double t1204;
  double t1228;
  double t1214;
  double t1219;
  double t1260;
  double t1258;
  double t1240;
  double t1244;
  double t1246;
  double t1247;
  double t1220;
  double t1227;
  double t1231;
  double t1234;
  double t1272;
  double t1253;
  double t1276;
  double t1280;
  double t1283;
  double t1294;
  double t1295;
  double t1296;
  double t1297;
  double t1298;
  double t1299;
  double t1300;
  double t1302;
  double t1208;
  double t1237;
  double t1238;
  double t1250;
  double t1251;
  double t1259;
  double t1267;
  double t1268;
  double t1315;
  double t1334;
  double t1335;
  double t1340;
  double t1356;
  double t1362;
  double t1363;
  double t1364;
  double t1384;
  double t1385;
  double t1386;
  double t1378;
  double t1379;
  double t1392;
  double t1393;
  double t1394;
  double t1395;
  double t1396;
  double t1399;
  double t1401;
  double t1404;
  double t1408;
  double t1409;
  double t1389;
  double t1390;
  double t1382;
  double t1387;
  double t1426;
  double t1427;
  double t1428;
  double t1423;
  double t1424;
  double t1425;
  double t1429;
  double t1430;
  double t1432;
  double t1433;
  double t1434;
  double t1435;
  double t1436;
  double t1437;
  double t1438;
  double t1357;
  double t1358;
  double t1360;
  double t1361;
  double t1366;
  double t1367;
  double t1368;
  double t1371;
  double t1344;
  double t1345;
  double t1372;
  double t1373;
  double t1374;
  double t1388;
  double t1410;
  double t1414;
  double t1468;
  double t1469;
  double t1470;
  double t1464;
  double t1465;
  double t1466;
  double t1416;
  t1204 = Cos(var1[4]);
  t1228 = Sin(var1[4]);
  t1214 = -1.*t1204;
  t1219 = 1. + t1214;
  t1260 = Cos(var1[3]);
  t1258 = Sin(var1[3]);
  t1240 = -0.0265*t1219;
  t1244 = -0.025229*t1204;
  t1246 = 0.07701400000000003*t1228;
  t1247 = t1240 + t1244 + t1246;
  t1220 = -0.2375*t1219;
  t1227 = -0.314514*t1204;
  t1231 = 0.0012709999999999978*t1228;
  t1234 = t1220 + t1227 + t1231;
  t1272 = Cos(var1[2]);
  t1253 = Sin(var1[2]);
  t1276 = t1260*t1204;
  t1280 = -1.*t1258*t1228;
  t1283 = t1276 + t1280;
  t1294 = 0.0265*t1204;
  t1295 = t1204*t1247;
  t1296 = 0.0695*t1228;
  t1297 = t1234*t1228;
  t1298 = t1294 + t1295 + t1296 + t1297;
  t1299 = -1.*t1204*t1258;
  t1300 = -1.*t1260*t1228;
  t1302 = t1299 + t1300;
  t1208 = -0.0695*t1204;
  t1237 = -1.*t1204*t1234;
  t1238 = 0.0265*t1228;
  t1250 = t1247*t1228;
  t1251 = t1208 + t1237 + t1238 + t1250;
  t1259 = t1204*t1258;
  t1267 = t1260*t1228;
  t1268 = t1259 + t1267;
  t1315 = -1.*t1253*t1283;
  t1334 = -1.*t1253*t1302;
  t1335 = t1272*t1283;
  t1340 = t1334 + t1335;
  t1356 = t1272*t1302;
  t1362 = -1.*t1260*t1204;
  t1363 = t1258*t1228;
  t1364 = t1362 + t1363;
  t1384 = -0.0695*t1258;
  t1385 = -1.*t1258*t1234;
  t1386 = t1260*t1247;
  t1378 = -1.*t1260;
  t1379 = 1. + t1378;
  t1392 = -0.0695*t1260;
  t1393 = -0.0265*t1258;
  t1394 = -1.*t1260*t1234;
  t1395 = -1.*t1258*t1247;
  t1396 = t1392 + t1393 + t1394 + t1395;
  t1399 = -0.0695*t1379;
  t1401 = 0.0265*t1258;
  t1404 = t1260*t1234;
  t1408 = t1258*t1247;
  t1409 = t1399 + t1401 + t1404 + t1408;
  t1389 = 0.0265*t1260;
  t1390 = t1389 + t1384 + t1385 + t1386;
  t1382 = -0.0265*t1379;
  t1387 = t1382 + t1384 + t1385 + t1386;
  t1426 = 0.07701400000000003*t1204;
  t1427 = -0.0012709999999999978*t1228;
  t1428 = t1426 + t1427;
  t1423 = 0.0012709999999999978*t1204;
  t1424 = t1423 + t1246;
  t1425 = -1.*t1204*t1424;
  t1429 = t1428*t1228;
  t1430 = t1294 + t1295 + t1425 + t1296 + t1429 + t1297;
  t1432 = 0.0695*t1204;
  t1433 = t1204*t1428;
  t1434 = t1204*t1234;
  t1435 = -0.0265*t1228;
  t1436 = -1.*t1247*t1228;
  t1437 = t1424*t1228;
  t1438 = t1432 + t1433 + t1434 + t1435 + t1436 + t1437;
  t1357 = t1253*t1283;
  t1358 = t1356 + t1357;
  t1360 = 0.19964*t1251*t1358;
  t1361 = t1253*t1302;
  t1366 = t1272*t1364;
  t1367 = t1361 + t1366;
  t1368 = 0.19964*t1298*t1367;
  t1371 = 0.19964*t1251*t1340;
  t1344 = t1272*t1268;
  t1345 = t1344 + t1315;
  t1372 = -1.*t1253*t1364;
  t1373 = t1356 + t1372;
  t1374 = 0.19964*t1298*t1373;
  t1388 = -1.*t1387*t1302;
  t1410 = -1.*t1409*t1283;
  t1414 = t1409*t1302;
  t1468 = t1260*t1428;
  t1469 = -1.*t1258*t1424;
  t1470 = t1468 + t1469;
  t1464 = t1258*t1428;
  t1465 = t1260*t1424;
  t1466 = t1464 + t1465;
  t1416 = t1387*t1364;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.19964*t1298*t1340 + 0.19964*t1251*t1345)*var2[0] + 0.5*(0.19964*t1251*(-1.*t1253*t1268 - 1.*t1272*t1283) + 0.19964*t1298*(-1.*t1272*t1302 + t1315))*var2[1])*var2[3];
  p_output1[3]=(0.5*(t1360 + t1368)*var2[0] + 0.5*(t1371 + t1374)*var2[1] + 0.5*(0.19964*t1298*(t1388 - 1.*t1268*t1390 - 1.*t1283*t1396 + t1410) + 0.19964*t1251*(t1283*t1390 + t1302*t1396 + t1414 + t1416))*var2[2])*var2[3];
  p_output1[4]=var2[3]*(0.5*(t1360 + t1368 + 0.19964*(t1253*t1268 + t1335)*t1430 + 0.19964*t1358*t1438)*var2[0] + 0.5*(t1371 + t1374 + 0.19964*t1345*t1430 + 0.19964*t1340*t1438)*var2[1] + 0.5*(0.19964*(t1302*t1387 + t1283*t1409)*t1430 + 0.19964*(-1.*t1283*t1387 - 1.*t1268*t1409)*t1438 + 0.19964*t1298*(t1388 + t1410 - 1.*t1268*t1466 - 1.*t1283*t1470) + 0.19964*t1251*(t1414 + t1416 + t1283*t1466 + t1302*t1470))*var2[2] + 0.5*(0.39928*t1251*t1430 + 0.39928*t1298*t1438)*var2[3] + 0.5*(0.015375074960000006*t1430 + 0.0002537424399999996*t1438)*var2[4]);
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

#include "Ce3_vec_L3_J4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L3_J4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
