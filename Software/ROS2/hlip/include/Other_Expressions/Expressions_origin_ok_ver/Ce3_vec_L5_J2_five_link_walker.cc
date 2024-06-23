/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:41:05 GMT-05:00
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
  double t1254;
  double t1224;
  double t1247;
  double t1255;
  double t1273;
  double t1196;
  double t1275;
  double t1279;
  double t1286;
  double t1253;
  double t1265;
  double t1267;
  double t1270;
  double t1294;
  double t1302;
  double t1341;
  double t1342;
  double t1344;
  double t1331;
  double t1333;
  double t1335;
  double t1336;
  double t1337;
  double t1338;
  double t1313;
  double t1314;
  double t1316;
  double t1324;
  double t1325;
  double t1326;
  double t1366;
  double t1363;
  double t1367;
  double t1358;
  double t1359;
  double t1400;
  double t1401;
  double t1404;
  double t1396;
  double t1397;
  double t1398;
  double t1349;
  double t1350;
  double t1351;
  double t1339;
  double t1345;
  double t1346;
  double t1422;
  double t1423;
  double t1424;
  double t1425;
  double t1426;
  double t1375;
  double t1376;
  double t1384;
  double t1385;
  double t1390;
  double t1411;
  double t1412;
  double t1413;
  double t1399;
  double t1405;
  double t1406;
  double t1454;
  double t1455;
  double t1456;
  double t1421;
  double t1427;
  double t1429;
  double t1430;
  double t1431;
  double t1475;
  double t1476;
  double t1477;
  double t1472;
  double t1473;
  double t1433;
  double t1436;
  double t1437;
  double t1438;
  double t1439;
  double t1440;
  double t1441;
  double t1442;
  double t1443;
  double t1445;
  double t1446;
  double t1449;
  double t1450;
  double t1451;
  double t1452;
  double t1458;
  double t1461;
  double t1463;
  double t1493;
  double t1494;
  double t1495;
  double t1489;
  double t1490;
  double t1491;
  double t1465;
  t1254 = Cos(var1[5]);
  t1224 = Cos(var1[6]);
  t1247 = Sin(var1[5]);
  t1255 = Sin(var1[6]);
  t1273 = Cos(var1[2]);
  t1196 = Sin(var1[2]);
  t1275 = t1254*t1224;
  t1279 = -1.*t1247*t1255;
  t1286 = t1275 + t1279;
  t1253 = -1.*t1224*t1247;
  t1265 = -1.*t1254*t1255;
  t1267 = t1253 + t1265;
  t1270 = -1.*t1196*t1267;
  t1294 = -1.*t1273*t1286;
  t1302 = t1270 + t1294;
  t1341 = -0.022659*t1224;
  t1342 = -0.007367999999999986*t1255;
  t1344 = t1341 + t1342;
  t1331 = -1.*t1224;
  t1333 = 1. + t1331;
  t1335 = -0.16*t1333;
  t1336 = -0.167368*t1224;
  t1337 = 0.022659*t1255;
  t1338 = t1335 + t1336 + t1337;
  t1313 = t1224*t1247;
  t1314 = t1254*t1255;
  t1316 = t1313 + t1314;
  t1324 = -1.*t1273*t1316;
  t1325 = -1.*t1196*t1286;
  t1326 = t1324 + t1325;
  t1366 = t1273*t1286;
  t1363 = -1.*t1196*t1316;
  t1367 = t1363 + t1366;
  t1358 = t1273*t1267;
  t1359 = t1358 + t1325;
  t1400 = -1.*t1247*t1344;
  t1401 = t1254*t1338;
  t1404 = t1400 + t1401;
  t1396 = t1254*t1344;
  t1397 = t1247*t1338;
  t1398 = t1396 + t1397;
  t1349 = -1.*t1224*t1344;
  t1350 = t1338*t1255;
  t1351 = t1349 + t1350;
  t1339 = t1224*t1338;
  t1345 = t1344*t1255;
  t1346 = t1339 + t1345;
  t1422 = -1.*t1254*t1224;
  t1423 = t1247*t1255;
  t1424 = t1422 + t1423;
  t1425 = t1273*t1424;
  t1426 = t1270 + t1425;
  t1375 = t1196*t1267;
  t1376 = t1375 + t1366;
  t1384 = t1273*t1316;
  t1385 = t1196*t1286;
  t1390 = t1384 + t1385;
  t1411 = -1.*t1404*t1267;
  t1412 = -1.*t1398*t1286;
  t1413 = t1411 + t1412;
  t1399 = t1398*t1316;
  t1405 = t1404*t1286;
  t1406 = t1399 + t1405;
  t1454 = -1.*t1254*t1344;
  t1455 = -1.*t1247*t1338;
  t1456 = t1454 + t1455;
  t1421 = 0.0033974904599999994*t1359;
  t1427 = -0.0011047579199999977*t1426;
  t1429 = t1421 + t1427;
  t1430 = 0.5*var2[6]*t1429;
  t1431 = 0.14994*t1351*t1359;
  t1475 = 0.022659*t1224;
  t1476 = 0.007367999999999986*t1255;
  t1477 = t1475 + t1476;
  t1472 = -0.007367999999999986*t1224;
  t1473 = t1472 + t1337;
  t1433 = 0.14994*t1346*t1426;
  t1436 = 0.29988*t1367*t1359;
  t1437 = 0.29988*t1359*t1426;
  t1438 = t1436 + t1437;
  t1439 = 0.5*var2[1]*t1438;
  t1440 = 0.14994*t1376*t1367;
  t1441 = 0.14994*t1359*t1390;
  t1442 = 0.14994*t1376*t1426;
  t1443 = t1196*t1424;
  t1445 = t1358 + t1443;
  t1446 = 0.14994*t1359*t1445;
  t1449 = t1440 + t1441 + t1442 + t1446;
  t1450 = 0.5*var2[0]*t1449;
  t1451 = 0.14994*t1359*t1413;
  t1452 = t1404*t1267;
  t1458 = t1398*t1286;
  t1461 = 0.14994*t1406*t1426;
  t1463 = -1.*t1398*t1267;
  t1493 = t1254*t1477;
  t1494 = -1.*t1247*t1473;
  t1495 = t1493 + t1494;
  t1489 = t1247*t1477;
  t1490 = t1254*t1473;
  t1491 = t1489 + t1490;
  t1465 = -1.*t1404*t1424;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(0.5*(0.14994*Power(t1359,2) + 0.14994*Power(t1367,2) + 0.14994*t1302*t1376 + 0.14994*t1326*t1390)*var2[0] + 0.5*(0.29988*t1302*t1359 + 0.29988*t1326*t1367)*var2[1] + 0.5*(0.14994*t1302*t1406 + 0.14994*t1326*t1413)*var2[2] + 0.5*(0.14994*t1302*t1346 + 0.14994*t1326*t1351)*var2[5] + 0.5*(-0.0011047579199999977*t1302 + 0.0033974904599999994*t1326)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(t1430 + t1439 + t1450 + 0.5*(t1451 + 0.14994*t1359*(t1316*t1404 + t1452 + t1286*t1456 + t1458) + t1461 + 0.14994*t1367*(-1.*t1286*t1404 - 1.*t1267*t1456 + t1463 + t1465))*var2[2] + 0.5*(t1431 + t1433)*var2[5]);
  p_output1[6]=var2[1]*(t1430 + t1439 + t1450 + 0.5*(t1451 + t1461 + 0.14994*t1367*(t1463 + t1465 - 1.*t1286*t1491 - 1.*t1267*t1495) + 0.14994*t1359*(t1452 + t1458 + t1316*t1491 + t1286*t1495))*var2[2] + 0.5*(t1431 + t1433 + 0.14994*t1359*(-1.*t1255*t1338 + t1224*t1344 + t1255*t1473 + t1224*t1477) + 0.14994*t1367*(t1339 + t1345 - 1.*t1224*t1473 + t1255*t1477))*var2[5]);
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

#include "Ce3_vec_L5_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L5_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
