/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:41:14 GMT-05:00
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
  double t1330;
  double t1346;
  double t1390;
  double t1383;
  double t1358;
  double t1359;
  double t1361;
  double t1367;
  double t1373;
  double t1374;
  double t1336;
  double t1347;
  double t1350;
  double t1399;
  double t1381;
  double t1401;
  double t1406;
  double t1410;
  double t1417;
  double t1419;
  double t1420;
  double t1421;
  double t1423;
  double t1425;
  double t1353;
  double t1376;
  double t1378;
  double t1384;
  double t1391;
  double t1396;
  double t1436;
  double t1426;
  double t1437;
  double t1457;
  double t1413;
  double t1475;
  double t1476;
  double t1478;
  double t1493;
  double t1496;
  double t1497;
  double t1500;
  double t1501;
  double t1502;
  double t1504;
  double t1505;
  double t1506;
  double t1530;
  double t1531;
  double t1532;
  double t1527;
  double t1528;
  double t1529;
  double t1533;
  double t1534;
  double t1536;
  double t1537;
  double t1538;
  double t1539;
  double t1540;
  double t1397;
  double t1414;
  double t1474;
  double t1479;
  double t1480;
  double t1481;
  double t1484;
  double t1485;
  double t1486;
  double t1487;
  double t1488;
  double t1489;
  double t1498;
  double t1508;
  double t1514;
  double t1571;
  double t1572;
  double t1573;
  double t1567;
  double t1568;
  double t1569;
  double t1518;
  t1330 = Cos(var1[6]);
  t1346 = Sin(var1[6]);
  t1390 = Cos(var1[5]);
  t1383 = Sin(var1[5]);
  t1358 = -1.*t1330;
  t1359 = 1. + t1358;
  t1361 = -0.16*t1359;
  t1367 = -0.167368*t1330;
  t1373 = 0.022659*t1346;
  t1374 = t1361 + t1367 + t1373;
  t1336 = -0.022659*t1330;
  t1347 = -0.007367999999999986*t1346;
  t1350 = t1336 + t1347;
  t1399 = Cos(var1[2]);
  t1381 = Sin(var1[2]);
  t1401 = t1390*t1330;
  t1406 = -1.*t1383*t1346;
  t1410 = t1401 + t1406;
  t1417 = t1330*t1374;
  t1419 = t1350*t1346;
  t1420 = t1417 + t1419;
  t1421 = -1.*t1330*t1383;
  t1423 = -1.*t1390*t1346;
  t1425 = t1421 + t1423;
  t1353 = -1.*t1330*t1350;
  t1376 = t1374*t1346;
  t1378 = t1353 + t1376;
  t1384 = t1330*t1383;
  t1391 = t1390*t1346;
  t1396 = t1384 + t1391;
  t1436 = -1.*t1381*t1410;
  t1426 = t1399*t1425;
  t1437 = t1426 + t1436;
  t1457 = -1.*t1381*t1425;
  t1413 = t1399*t1410;
  t1475 = -1.*t1390*t1330;
  t1476 = t1383*t1346;
  t1478 = t1475 + t1476;
  t1493 = -1.*t1383*t1350;
  t1496 = t1390*t1374;
  t1497 = t1493 + t1496;
  t1500 = -1.*t1390*t1350;
  t1501 = -1.*t1383*t1374;
  t1502 = t1500 + t1501;
  t1504 = t1390*t1350;
  t1505 = t1383*t1374;
  t1506 = t1504 + t1505;
  t1530 = 0.022659*t1330;
  t1531 = 0.007367999999999986*t1346;
  t1532 = t1530 + t1531;
  t1527 = -0.007367999999999986*t1330;
  t1528 = t1527 + t1373;
  t1529 = -1.*t1330*t1528;
  t1533 = t1532*t1346;
  t1534 = t1417 + t1529 + t1419 + t1533;
  t1536 = t1330*t1350;
  t1537 = t1330*t1532;
  t1538 = -1.*t1374*t1346;
  t1539 = t1528*t1346;
  t1540 = t1536 + t1537 + t1538 + t1539;
  t1397 = -1.*t1381*t1396;
  t1414 = t1397 + t1413;
  t1474 = 0.14994*t1378*t1437;
  t1479 = t1399*t1478;
  t1480 = t1457 + t1479;
  t1481 = 0.14994*t1420*t1480;
  t1484 = t1381*t1425;
  t1485 = t1484 + t1413;
  t1486 = 0.14994*t1378*t1485;
  t1487 = t1381*t1478;
  t1488 = t1426 + t1487;
  t1489 = 0.14994*t1420*t1488;
  t1498 = t1497*t1425;
  t1508 = t1506*t1410;
  t1514 = -1.*t1506*t1425;
  t1571 = t1390*t1532;
  t1572 = -1.*t1383*t1528;
  t1573 = t1571 + t1572;
  t1567 = t1383*t1532;
  t1568 = t1390*t1528;
  t1569 = t1567 + t1568;
  t1518 = -1.*t1497*t1478;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.14994*t1378*t1414 + 0.14994*t1420*t1437)*var2[0] + 0.5*(0.14994*t1378*(-1.*t1396*t1399 + t1436) + 0.14994*t1420*(-1.*t1399*t1410 + t1457))*var2[1])*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(0.5*(t1486 + t1489)*var2[0] + 0.5*(t1474 + t1481)*var2[1] + 0.5*(0.14994*t1420*(t1396*t1497 + t1498 + t1410*t1502 + t1508) + 0.14994*t1378*(-1.*t1410*t1497 - 1.*t1425*t1502 + t1514 + t1518))*var2[2])*var2[5];
  p_output1[6]=var2[5]*(0.5*(t1486 + t1489 + 0.14994*(t1396*t1399 + t1381*t1410)*t1534 + 0.14994*t1485*t1540)*var2[0] + 0.5*(t1474 + t1481 + 0.14994*t1414*t1534 + 0.14994*t1437*t1540)*var2[1] + 0.5*(0.14994*(-1.*t1425*t1497 - 1.*t1410*t1506)*t1534 + 0.14994*(t1410*t1497 + t1396*t1506)*t1540 + 0.14994*t1420*(t1498 + t1508 + t1396*t1569 + t1410*t1573) + 0.14994*t1378*(t1514 + t1518 - 1.*t1410*t1569 - 1.*t1425*t1573))*var2[2] + 0.5*(0.29988*t1378*t1534 + 0.29988*t1420*t1540)*var2[5] + 0.5*(0.0033974904599999994*t1534 - 0.0011047579199999977*t1540)*var2[6]);
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

#include "Ce3_vec_L5_J6_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L5_J6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
