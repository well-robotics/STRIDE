/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:34:01 GMT-05:00
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
  double t1455;
  double t1437;
  double t1438;
  double t1467;
  double t1482;
  double t1434;
  double t1483;
  double t1488;
  double t1489;
  double t1510;
  double t1511;
  double t1448;
  double t1472;
  double t1474;
  double t1477;
  double t1491;
  double t1494;
  double t1518;
  double t1519;
  double t1520;
  double t1521;
  double t1512;
  double t1513;
  double t1514;
  double t1515;
  double t1497;
  double t1499;
  double t1500;
  double t1501;
  double t1502;
  double t1504;
  double t1533;
  double t1534;
  double t1536;
  double t1537;
  double t1538;
  double t1555;
  double t1556;
  double t1563;
  double t1564;
  double t1565;
  double t1566;
  double t1567;
  double t1557;
  double t1558;
  double t1559;
  double t1560;
  double t1561;
  double t1509;
  double t1516;
  double t1517;
  double t1522;
  double t1523;
  double t1525;
  double t1526;
  double t1527;
  double t1528;
  double t1529;
  double t1580;
  double t1581;
  double t1582;
  double t1583;
  double t1584;
  double t1542;
  double t1571;
  double t1572;
  double t1573;
  double t1562;
  double t1568;
  double t1569;
  double t1609;
  double t1610;
  double t1611;
  double t1612;
  double t1613;
  double t1606;
  double t1607;
  double t1579;
  double t1585;
  double t1586;
  double t1587;
  double t1588;
  double t1637;
  double t1638;
  double t1631;
  double t1632;
  double t1633;
  double t1589;
  double t1592;
  double t1593;
  double t1594;
  double t1595;
  double t1596;
  double t1597;
  double t1598;
  double t1599;
  double t1600;
  double t1601;
  double t1602;
  double t1603;
  double t1604;
  double t1605;
  double t1615;
  double t1618;
  double t1620;
  double t1652;
  double t1653;
  double t1654;
  double t1648;
  double t1649;
  double t1650;
  double t1622;
  t1455 = Cos(var1[5]);
  t1437 = Cos(var1[6]);
  t1438 = Sin(var1[5]);
  t1467 = Sin(var1[6]);
  t1482 = Cos(var1[2]);
  t1434 = Sin(var1[2]);
  t1483 = t1455*t1437;
  t1488 = -1.*t1438*t1467;
  t1489 = t1483 + t1488;
  t1510 = -1.*t1437;
  t1511 = 1. + t1510;
  t1448 = t1437*t1438;
  t1472 = t1455*t1467;
  t1474 = t1448 + t1472;
  t1477 = -1.*t1434*t1474;
  t1491 = t1482*t1489;
  t1494 = t1477 + t1491;
  t1518 = -0.2375*t1511;
  t1519 = -0.314506*t1437;
  t1520 = -0.0012740000000000008*t1467;
  t1521 = t1518 + t1519 + t1520;
  t1512 = -0.0265*t1511;
  t1513 = -0.025226*t1437;
  t1514 = -0.07700600000000002*t1467;
  t1515 = t1512 + t1513 + t1514;
  t1497 = -1.*t1437*t1438;
  t1499 = -1.*t1455*t1467;
  t1500 = t1497 + t1499;
  t1501 = t1482*t1500;
  t1502 = -1.*t1434*t1489;
  t1504 = t1501 + t1502;
  t1533 = t1434*t1500;
  t1534 = t1533 + t1491;
  t1536 = t1482*t1474;
  t1537 = t1434*t1489;
  t1538 = t1536 + t1537;
  t1555 = -1.*t1455;
  t1556 = 1. + t1555;
  t1563 = -0.0695*t1556;
  t1564 = -0.0265*t1438;
  t1565 = -1.*t1438*t1515;
  t1566 = t1455*t1521;
  t1567 = t1563 + t1564 + t1565 + t1566;
  t1557 = -0.0265*t1556;
  t1558 = 0.0695*t1438;
  t1559 = t1455*t1515;
  t1560 = t1438*t1521;
  t1561 = t1557 + t1558 + t1559 + t1560;
  t1509 = -0.0265*t1437;
  t1516 = -1.*t1437*t1515;
  t1517 = 0.0695*t1467;
  t1522 = t1521*t1467;
  t1523 = t1509 + t1516 + t1517 + t1522;
  t1525 = 0.0695*t1437;
  t1526 = t1437*t1521;
  t1527 = 0.0265*t1467;
  t1528 = t1515*t1467;
  t1529 = t1525 + t1526 + t1527 + t1528;
  t1580 = -1.*t1455*t1437;
  t1581 = t1438*t1467;
  t1582 = t1580 + t1581;
  t1583 = t1434*t1582;
  t1584 = t1501 + t1583;
  t1542 = -1.*t1434*t1500;
  t1571 = -1.*t1567*t1500;
  t1572 = -1.*t1561*t1489;
  t1573 = t1571 + t1572;
  t1562 = t1561*t1474;
  t1568 = t1567*t1489;
  t1569 = t1562 + t1568;
  t1609 = -0.0265*t1455;
  t1610 = -0.0695*t1438;
  t1611 = -1.*t1455*t1515;
  t1612 = -1.*t1438*t1521;
  t1613 = t1609 + t1610 + t1611 + t1612;
  t1606 = 0.0695*t1455;
  t1607 = t1606 + t1564 + t1565 + t1566;
  t1579 = -0.0002543413600000002*t1534;
  t1585 = -0.015373477840000005*t1584;
  t1586 = t1579 + t1585;
  t1587 = 0.5*var2[6]*t1586;
  t1588 = 0.19964*t1523*t1534;
  t1637 = -0.07700600000000002*t1437;
  t1638 = t1637 + t1520;
  t1631 = -0.0012740000000000008*t1437;
  t1632 = 0.07700600000000002*t1467;
  t1633 = t1631 + t1632;
  t1589 = 0.19964*t1529*t1584;
  t1592 = 0.39928*t1534*t1538;
  t1593 = 0.39928*t1534*t1584;
  t1594 = t1592 + t1593;
  t1595 = 0.5*var2[0]*t1594;
  t1596 = 0.19964*t1534*t1494;
  t1597 = 0.19964*t1504*t1538;
  t1598 = t1482*t1582;
  t1599 = t1542 + t1598;
  t1600 = 0.19964*t1534*t1599;
  t1601 = 0.19964*t1504*t1584;
  t1602 = t1596 + t1597 + t1600 + t1601;
  t1603 = 0.5*var2[1]*t1602;
  t1604 = 0.19964*t1534*t1573;
  t1605 = t1567*t1500;
  t1615 = t1561*t1489;
  t1618 = 0.19964*t1569*t1584;
  t1620 = -1.*t1561*t1500;
  t1652 = -1.*t1438*t1638;
  t1653 = t1455*t1633;
  t1654 = t1652 + t1653;
  t1648 = t1455*t1638;
  t1649 = t1438*t1633;
  t1650 = t1648 + t1649;
  t1622 = -1.*t1567*t1582;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(0.5*(0.39928*t1504*t1534 + 0.39928*t1494*t1538)*var2[0] + 0.5*(0.19964*Power(t1494,2) + 0.19964*Power(t1504,2) + 0.19964*(-1.*t1474*t1482 + t1502)*t1538 + 0.19964*t1534*(-1.*t1482*t1489 + t1542))*var2[1] + 0.5*(0.19964*t1504*t1569 + 0.19964*t1494*t1573)*var2[2] + 0.5*(0.19964*t1494*t1523 + 0.19964*t1504*t1529)*var2[5] + 0.5*(-0.0002543413600000002*t1494 - 0.015373477840000005*t1504)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(t1587 + t1595 + t1603 + 0.5*(t1604 + 0.19964*t1534*(t1605 + t1474*t1607 + t1489*t1613 + t1615) + t1618 + 0.19964*t1538*(-1.*t1489*t1607 - 1.*t1500*t1613 + t1620 + t1622))*var2[2] + 0.5*(t1588 + t1589)*var2[5]);
  p_output1[6]=var2[0]*(t1587 + t1595 + t1603 + 0.5*(t1604 + t1618 + 0.19964*t1534*(t1605 + t1615 + t1474*t1650 + t1489*t1654) + 0.19964*t1538*(t1620 + t1622 - 1.*t1489*t1650 - 1.*t1500*t1654))*var2[2] + 0.5*(t1588 + t1589 + 0.19964*t1538*(t1525 + t1526 + t1527 + t1528 + t1467*t1633 - 1.*t1437*t1638) + 0.19964*t1534*(0.0265*t1437 - 0.0695*t1467 + t1437*t1515 - 1.*t1467*t1521 + t1437*t1633 + t1467*t1638))*var2[5]);
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

#include "Ce3_vec_L5_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L5_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
