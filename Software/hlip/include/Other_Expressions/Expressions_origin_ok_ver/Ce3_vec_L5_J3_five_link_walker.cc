/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:41:08 GMT-05:00
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
  double t1314;
  double t1302;
  double t1305;
  double t1326;
  double t1347;
  double t1350;
  double t1351;
  double t1354;
  double t1357;
  double t1358;
  double t1359;
  double t1361;
  double t1366;
  double t1336;
  double t1337;
  double t1341;
  double t1335;
  double t1373;
  double t1374;
  double t1375;
  double t1275;
  double t1377;
  double t1378;
  double t1379;
  double t1313;
  double t1328;
  double t1329;
  double t1353;
  double t1367;
  double t1370;
  double t1376;
  double t1381;
  double t1383;
  double t1342;
  double t1392;
  double t1393;
  double t1396;
  double t1427;
  double t1434;
  double t1435;
  double t1425;
  double t1426;
  double t1436;
  double t1437;
  double t1440;
  double t1442;
  double t1443;
  double t1445;
  double t1446;
  double t1451;
  double t1452;
  double t1453;
  double t1457;
  double t1330;
  double t1346;
  double t1405;
  double t1385;
  double t1390;
  double t1391;
  double t1483;
  double t1484;
  double t1505;
  double t1506;
  double t1507;
  double t1509;
  double t1510;
  double t1514;
  double t1515;
  double t1516;
  double t1508;
  double t1511;
  double t1512;
  double t1461;
  double t1462;
  double t1464;
  double t1513;
  double t1517;
  double t1518;
  double t1467;
  double t1468;
  double t1469;
  double t1520;
  double t1521;
  double t1522;
  double t1474;
  double t1476;
  double t1478;
  double t1479;
  double t1485;
  double t1487;
  double t1488;
  double t1489;
  double t1490;
  double t1492;
  double t1493;
  t1314 = Cos(var1[5]);
  t1302 = Cos(var1[6]);
  t1305 = Sin(var1[5]);
  t1326 = Sin(var1[6]);
  t1347 = -0.022659*t1302;
  t1350 = -0.007367999999999986*t1326;
  t1351 = t1347 + t1350;
  t1354 = -1.*t1302;
  t1357 = 1. + t1354;
  t1358 = -0.16*t1357;
  t1359 = -0.167368*t1302;
  t1361 = 0.022659*t1326;
  t1366 = t1358 + t1359 + t1361;
  t1336 = t1314*t1302;
  t1337 = -1.*t1305*t1326;
  t1341 = t1336 + t1337;
  t1335 = Sin(var1[2]);
  t1373 = t1302*t1305;
  t1374 = t1314*t1326;
  t1375 = t1373 + t1374;
  t1275 = Cos(var1[2]);
  t1377 = -1.*t1305*t1351;
  t1378 = t1314*t1366;
  t1379 = t1377 + t1378;
  t1313 = -1.*t1302*t1305;
  t1328 = -1.*t1314*t1326;
  t1329 = t1313 + t1328;
  t1353 = t1314*t1351;
  t1367 = t1305*t1366;
  t1370 = t1353 + t1367;
  t1376 = t1370*t1375;
  t1381 = t1379*t1341;
  t1383 = t1376 + t1381;
  t1342 = -1.*t1335*t1341;
  t1392 = -1.*t1379*t1329;
  t1393 = -1.*t1370*t1341;
  t1396 = t1392 + t1393;
  t1427 = -1.*t1314*t1351;
  t1434 = -1.*t1305*t1366;
  t1435 = t1427 + t1434;
  t1425 = t1379*t1329;
  t1426 = t1379*t1375;
  t1436 = t1435*t1341;
  t1437 = t1370*t1341;
  t1440 = t1425 + t1426 + t1436 + t1437;
  t1442 = -1.*t1435*t1329;
  t1443 = -1.*t1370*t1329;
  t1445 = -1.*t1379*t1341;
  t1446 = -1.*t1314*t1302;
  t1451 = t1305*t1326;
  t1452 = t1446 + t1451;
  t1453 = -1.*t1379*t1452;
  t1457 = t1442 + t1443 + t1445 + t1453;
  t1330 = t1275*t1329;
  t1346 = t1330 + t1342;
  t1405 = -1.*t1335*t1329;
  t1385 = -1.*t1335*t1375;
  t1390 = t1275*t1341;
  t1391 = t1385 + t1390;
  t1483 = t1335*t1329;
  t1484 = t1483 + t1390;
  t1505 = 0.022659*t1302;
  t1506 = 0.007367999999999986*t1326;
  t1507 = t1505 + t1506;
  t1509 = -0.007367999999999986*t1302;
  t1510 = t1509 + t1361;
  t1514 = t1314*t1507;
  t1515 = -1.*t1305*t1510;
  t1516 = t1514 + t1515;
  t1508 = t1305*t1507;
  t1511 = t1314*t1510;
  t1512 = t1508 + t1511;
  t1461 = t1302*t1366;
  t1462 = t1351*t1326;
  t1464 = t1461 + t1462;
  t1513 = t1512*t1375;
  t1517 = t1516*t1341;
  t1518 = t1425 + t1513 + t1437 + t1517;
  t1467 = -1.*t1302*t1351;
  t1468 = t1366*t1326;
  t1469 = t1467 + t1468;
  t1520 = -1.*t1516*t1329;
  t1521 = -1.*t1512*t1341;
  t1522 = t1443 + t1520 + t1521 + t1453;
  t1474 = 0.14994*t1346*t1396;
  t1476 = t1275*t1452;
  t1478 = t1405 + t1476;
  t1479 = 0.14994*t1383*t1478;
  t1485 = 0.14994*t1484*t1396;
  t1487 = t1335*t1452;
  t1488 = t1330 + t1487;
  t1489 = 0.14994*t1383*t1488;
  t1490 = t1275*t1375;
  t1492 = t1335*t1341;
  t1493 = t1490 + t1492;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.14994*t1346*t1383 + 0.14994*t1391*t1396)*var2[0] + 0.5*(0.14994*(t1342 - 1.*t1275*t1375)*t1396 + 0.14994*t1383*(-1.*t1275*t1341 + t1405))*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[2]*(0.5*(0.14994*t1440*t1484 + t1485 + t1489 + 0.14994*t1457*t1493)*var2[0] + 0.5*(0.14994*t1346*t1440 + 0.14994*t1391*t1457 + t1474 + t1479)*var2[1] + 0.5*(0.29988*t1383*t1440 + 0.29988*t1396*t1457)*var2[2] + 0.5*(0.14994*t1440*t1464 + 0.14994*t1457*t1469)*var2[5] + 0.5*(-0.0011047579199999977*t1440 + 0.0033974904599999994*t1457)*var2[6]);
  p_output1[6]=var2[2]*(0.5*(t1485 + t1489 + 0.14994*t1484*t1518 + 0.14994*t1493*t1522)*var2[0] + 0.5*(t1474 + t1479 + 0.14994*t1346*t1518 + 0.14994*t1391*t1522)*var2[1] + 0.5*(0.29988*t1383*t1518 + 0.29988*t1396*t1522)*var2[2] + 0.5*(0.14994*t1396*(t1461 + t1462 + t1326*t1507 - 1.*t1302*t1510) + 0.14994*t1383*(t1302*t1351 - 1.*t1326*t1366 + t1302*t1507 + t1326*t1510) + 0.14994*t1464*t1518 + 0.14994*t1469*t1522)*var2[5] + 0.5*(-0.0011047579199999977*t1518 + 0.0033974904599999994*t1522)*var2[6]);
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

#include "Ce3_vec_L5_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L5_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
