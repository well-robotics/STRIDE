/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:35 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t1388;
  double t1385;
  double t1386;
  double t1389;
  double t1437;
  double t1340;
  double t1387;
  double t1402;
  double t1435;
  double t1458;
  double t1464;
  double t1465;
  double t1466;
  double t1467;
  double t1436;
  double t1438;
  double t1444;
  double t1445;
  double t1446;
  double t1447;
  double t1449;
  double t1455;
  double t1456;
  double t1485;
  double t1486;
  double t1487;
  double t1471;
  double t1472;
  double t1473;
  double t1474;
  double t1475;
  double t1476;
  double t1489;
  double t1490;
  double t1491;
  double t1493;
  double t1494;
  double t1495;
  double t1496;
  double t1497;
  double t1498;
  double t1515;
  double t1516;
  double t1531;
  double t1532;
  double t1533;
  double t1535;
  double t1536;
  double t1537;
  double t1541;
  double t1542;
  double t1543;
  double t1518;
  double t1519;
  double t1520;
  double t1508;
  double t1509;
  double t1510;
  double t1468;
  double t1469;
  double t1470;
  double t1478;
  double t1479;
  double t1480;
  double t1481;
  double t1488;
  double t1512;
  double t1513;
  double t1503;
  double t1504;
  double t1505;
  double t1506;
  double t1507;
  double t1511;
  double t1514;
  double t1517;
  double t1521;
  double t1522;
  double t1523;
  double t1525;
  double t1526;
  double t1527;
  double t1528;
  double t1529;
  double t1534;
  double t1538;
  double t1539;
  double t1544;
  double t1545;
  double t1546;
  double t1547;
  double t1548;
  double t1550;
  double t1551;
  double t1552;
  double t1554;
  double t1555;
  double t1556;
  double t1557;
  double t1558;
  double t1576;
  double t1577;
  double t1578;
  double t1579;
  double t1580;
  double t1581;
  double t1530;
  double t1540;
  double t1549;
  double t1553;
  double t1559;
  double t1560;
  double t1565;
  double t1566;
  double t1567;
  double t1568;
  double t1457;
  double t1477;
  double t1482;
  double t1483;
  double t1586;
  double t1587;
  double t1588;
  t1388 = Cos(var1[5]);
  t1385 = Cos(var1[6]);
  t1386 = Sin(var1[5]);
  t1389 = Sin(var1[6]);
  t1437 = Sin(var1[2]);
  t1340 = Cos(var1[2]);
  t1387 = -1.*t1385*t1386;
  t1402 = -1.*t1388*t1389;
  t1435 = t1387 + t1402;
  t1458 = -1.*t1385;
  t1464 = 1. + t1458;
  t1465 = 0.5*t1464;
  t1466 = 0.671885*t1385;
  t1467 = t1465 + t1466;
  t1436 = t1340*t1435;
  t1438 = -1.*t1388*t1385;
  t1444 = t1386*t1389;
  t1445 = t1438 + t1444;
  t1446 = t1437*t1445;
  t1447 = t1436 + t1446;
  t1449 = -1.*t1388*t1437;
  t1455 = -1.*t1340*t1386;
  t1456 = t1449 + t1455;
  t1485 = t1340*t1388;
  t1486 = -1.*t1437*t1386;
  t1487 = t1485 + t1486;
  t1471 = t1437*t1435;
  t1472 = t1388*t1385;
  t1473 = -1.*t1386*t1389;
  t1474 = t1472 + t1473;
  t1475 = t1340*t1474;
  t1476 = t1471 + t1475;
  t1489 = t1388*t1437;
  t1490 = t1340*t1386;
  t1491 = t1489 + t1490;
  t1493 = t1385*t1386;
  t1494 = t1388*t1389;
  t1495 = t1493 + t1494;
  t1496 = t1340*t1495;
  t1497 = t1437*t1474;
  t1498 = t1496 + t1497;
  t1515 = -1.*t1437*t1474;
  t1516 = t1436 + t1515;
  t1531 = t1467*t1386;
  t1532 = 0.171885*t1388*t1389;
  t1533 = t1531 + t1532;
  t1535 = t1388*t1467;
  t1536 = -0.171885*t1386*t1389;
  t1537 = t1535 + t1536;
  t1541 = -1.*t1467*t1386;
  t1542 = -0.171885*t1388*t1389;
  t1543 = t1541 + t1542;
  t1518 = -1.*t1437*t1435;
  t1519 = t1340*t1445;
  t1520 = t1518 + t1519;
  t1508 = -1.*t1340*t1388;
  t1509 = t1437*t1386;
  t1510 = t1508 + t1509;
  t1468 = t1467*t1389;
  t1469 = -0.171885*t1385*t1389;
  t1470 = t1468 + t1469;
  t1478 = t1467*t1385;
  t1479 = Power(t1389,2);
  t1480 = 0.171885*t1479;
  t1481 = t1478 + t1480;
  t1488 = 6.8522*t1456*t1487;
  t1512 = -1.*t1437*t1495;
  t1513 = t1512 + t1475;
  t1503 = Power(t1456,2);
  t1504 = 3.4261*t1503;
  t1505 = 3.4261*t1456*t1491;
  t1506 = Power(t1487,2);
  t1507 = 3.4261*t1506;
  t1511 = 3.4261*t1487*t1510;
  t1514 = 0.85216*t1476*t1513;
  t1517 = 0.85216*t1516*t1498;
  t1521 = 0.85216*t1476*t1520;
  t1522 = 0.85216*t1516*t1447;
  t1523 = t1504 + t1505 + t1507 + t1511 + t1514 + t1517 + t1521 + t1522;
  t1525 = Power(t1388,2);
  t1526 = 0.1494*t1525;
  t1527 = Power(t1386,2);
  t1528 = 0.1494*t1527;
  t1529 = t1526 + t1528;
  t1534 = -1.*t1533*t1474;
  t1538 = -1.*t1435*t1537;
  t1539 = t1534 + t1538;
  t1544 = t1543*t1474;
  t1545 = t1533*t1474;
  t1546 = t1435*t1537;
  t1547 = t1495*t1537;
  t1548 = t1544 + t1545 + t1546 + t1547;
  t1550 = t1533*t1495;
  t1551 = t1474*t1537;
  t1552 = t1550 + t1551;
  t1554 = -1.*t1435*t1543;
  t1555 = -1.*t1435*t1533;
  t1556 = -1.*t1474*t1537;
  t1557 = -1.*t1537*t1445;
  t1558 = t1554 + t1555 + t1556 + t1557;
  t1576 = 3.4261*t1510*t1529;
  t1577 = 0.85216*t1516*t1539;
  t1578 = 0.85216*t1516*t1548;
  t1579 = 0.85216*t1552*t1520;
  t1580 = 0.85216*t1513*t1558;
  t1581 = t1576 + t1577 + t1578 + t1579 + t1580;
  t1530 = 3.4261*t1456*t1529;
  t1540 = 0.85216*t1476*t1539;
  t1549 = 0.85216*t1476*t1548;
  t1553 = 0.85216*t1552*t1447;
  t1559 = 0.85216*t1498*t1558;
  t1560 = t1530 + t1540 + t1549 + t1553 + t1559;
  t1565 = 0.51185934*t1510;
  t1566 = 0.85216*t1470*t1516;
  t1567 = 0.85216*t1481*t1520;
  t1568 = t1565 + t1566 + t1567;
  t1457 = 0.51185934*t1456;
  t1477 = 0.85216*t1470*t1476;
  t1482 = 0.85216*t1481*t1447;
  t1483 = t1457 + t1477 + t1482;
  t1586 = 0.85216*t1481*t1548;
  t1587 = 0.85216*t1470*t1558;
  t1588 = t1586 + t1587;
  p_output1[0]=var2[5]*(-0.5*(1.70432*t1447*t1476 + t1488 + 6.8522*t1487*t1491 + 1.70432*t1476*t1498)*var2[0] - 0.5*t1523*var2[1] - 0.5*t1560*var2[2] - 0.5*t1483*var2[5] - 0.0732367608*t1447*var2[6]);
  p_output1[1]=var2[5]*(-0.5*t1523*var2[0] - 0.5*(t1488 + 6.8522*t1456*t1510 + 1.70432*t1513*t1516 + 1.70432*t1516*t1520)*var2[1] - 0.5*t1581*var2[2] - 0.5*t1568*var2[5] - 0.0732367608*t1520*var2[6]);
  p_output1[2]=var2[5]*(-0.5*t1560*var2[0] - 0.5*t1581*var2[1] - 0.5*(1.70432*t1548*t1552 + 1.70432*t1539*t1558)*var2[2] - 0.5*t1588*var2[5] - 0.0732367608*t1548*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t1483*var2[0] - 0.5*t1568*var2[1] - 0.5*t1588*var2[2])*var2[5];
  p_output1[6]=(-0.0732367608*t1447*var2[0] - 0.0732367608*t1520*var2[1] - 0.0732367608*t1548*var2[2])*var2[5];
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

#include "Ce2_vec6_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
