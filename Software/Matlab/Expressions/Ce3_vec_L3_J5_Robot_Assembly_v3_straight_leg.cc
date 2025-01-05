/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:37 GMT-05:00
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
  double t1238;
  double t1220;
  double t1227;
  double t1240;
  double t1251;
  double t1208;
  double t1259;
  double t1267;
  double t1271;
  double t1293;
  double t1294;
  double t1295;
  double t1231;
  double t1244;
  double t1246;
  double t1297;
  double t1326;
  double t1332;
  double t1335;
  double t1354;
  double t1363;
  double t1366;
  double t1367;
  double t1404;
  double t1408;
  double t1401;
  double t1411;
  double t1412;
  double t1413;
  double t1414;
  double t1415;
  double t1416;
  double t1417;
  double t1418;
  double t1419;
  double t1420;
  double t1393;
  double t1398;
  double t1428;
  double t1429;
  double t1430;
  double t1431;
  double t1432;
  double t1434;
  double t1435;
  double t1437;
  double t1438;
  double t1439;
  double t1423;
  double t1426;
  double t1399;
  double t1421;
  double t1459;
  double t1460;
  double t1461;
  double t1455;
  double t1456;
  double t1356;
  double t1357;
  double t1358;
  double t1361;
  double t1369;
  double t1370;
  double t1373;
  double t1375;
  double t1376;
  double t1382;
  double t1384;
  double t1386;
  double t1389;
  double t1391;
  double t1392;
  double t1422;
  double t1440;
  double t1444;
  double t1483;
  double t1484;
  double t1485;
  double t1479;
  double t1480;
  double t1481;
  double t1446;
  t1238 = Cos(var1[3]);
  t1220 = Cos(var1[4]);
  t1227 = Sin(var1[3]);
  t1240 = Sin(var1[4]);
  t1251 = Cos(var1[2]);
  t1208 = Sin(var1[2]);
  t1259 = t1238*t1220;
  t1267 = -1.*t1227*t1240;
  t1271 = t1259 + t1267;
  t1293 = -1.*t1220*t1227;
  t1294 = -1.*t1238*t1240;
  t1295 = t1293 + t1294;
  t1231 = t1220*t1227;
  t1244 = t1238*t1240;
  t1246 = t1231 + t1244;
  t1297 = -1.*t1208*t1271;
  t1326 = -1.*t1208*t1295;
  t1332 = t1251*t1271;
  t1335 = t1326 + t1332;
  t1354 = t1251*t1295;
  t1363 = -1.*t1238*t1220;
  t1366 = t1227*t1240;
  t1367 = t1363 + t1366;
  t1404 = -1.*t1220;
  t1408 = 1. + t1404;
  t1401 = -0.0695*t1227;
  t1411 = -0.2375*t1408;
  t1412 = -0.314514*t1220;
  t1413 = 0.0012709999999999978*t1240;
  t1414 = t1411 + t1412 + t1413;
  t1415 = -1.*t1227*t1414;
  t1416 = -0.0265*t1408;
  t1417 = -0.025229*t1220;
  t1418 = 0.07701400000000003*t1240;
  t1419 = t1416 + t1417 + t1418;
  t1420 = t1238*t1419;
  t1393 = -1.*t1238;
  t1398 = 1. + t1393;
  t1428 = -0.0695*t1238;
  t1429 = -0.0265*t1227;
  t1430 = -1.*t1238*t1414;
  t1431 = -1.*t1227*t1419;
  t1432 = t1428 + t1429 + t1430 + t1431;
  t1434 = -0.0695*t1398;
  t1435 = 0.0265*t1227;
  t1437 = t1238*t1414;
  t1438 = t1227*t1419;
  t1439 = t1434 + t1435 + t1437 + t1438;
  t1423 = 0.0265*t1238;
  t1426 = t1423 + t1401 + t1415 + t1420;
  t1399 = -0.0265*t1398;
  t1421 = t1399 + t1401 + t1415 + t1420;
  t1459 = 0.07701400000000003*t1220;
  t1460 = -0.0012709999999999978*t1240;
  t1461 = t1459 + t1460;
  t1455 = 0.0012709999999999978*t1220;
  t1456 = t1455 + t1418;
  t1356 = t1208*t1271;
  t1357 = t1354 + t1356;
  t1358 = 0.015375074960000006*t1357;
  t1361 = t1208*t1295;
  t1369 = t1251*t1367;
  t1370 = t1361 + t1369;
  t1373 = 0.0002537424399999996*t1370;
  t1375 = t1358 + t1373;
  t1376 = 0.5*var2[0]*t1375;
  t1382 = 0.015375074960000006*t1335;
  t1384 = -1.*t1208*t1367;
  t1386 = t1354 + t1384;
  t1389 = 0.0002537424399999996*t1386;
  t1391 = t1382 + t1389;
  t1392 = 0.5*var2[1]*t1391;
  t1422 = -1.*t1421*t1295;
  t1440 = -1.*t1439*t1271;
  t1444 = t1439*t1295;
  t1483 = t1238*t1461;
  t1484 = -1.*t1227*t1456;
  t1485 = t1483 + t1484;
  t1479 = t1227*t1461;
  t1480 = t1238*t1456;
  t1481 = t1479 + t1480;
  t1446 = t1421*t1367;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.015375074960000006*(t1246*t1251 + t1297) + 0.0002537424399999996*t1335)*var2[0] + 0.5*(0.015375074960000006*(-1.*t1208*t1246 - 1.*t1251*t1271) + 0.0002537424399999996*(-1.*t1251*t1295 + t1297))*var2[1])*var2[4];
  p_output1[3]=(t1376 + t1392 + 0.5*(0.0002537424399999996*(t1422 - 1.*t1246*t1426 - 1.*t1271*t1432 + t1440) + 0.015375074960000006*(t1271*t1426 + t1295*t1432 + t1444 + t1446))*var2[2])*var2[4];
  p_output1[4]=(t1376 + t1392 + 0.5*(0.0002537424399999996*(t1422 + t1440 - 1.*t1246*t1481 - 1.*t1271*t1485) + 0.015375074960000006*(t1444 + t1446 + t1271*t1481 + t1295*t1485))*var2[2] + 0.5*(0.0002537424399999996*(0.0695*t1220 - 0.0265*t1240 + t1220*t1414 - 1.*t1240*t1419 + t1240*t1456 + t1220*t1461) + 0.015375074960000006*(0.0265*t1220 + 0.0695*t1240 + t1240*t1414 + t1220*t1419 - 1.*t1220*t1456 + t1240*t1461))*var2[3])*var2[4];
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

#include "Ce3_vec_L3_J5_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L3_J5_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
