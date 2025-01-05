/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:34:04 GMT-05:00
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
  double t1494;
  double t1472;
  double t1483;
  double t1495;
  double t1507;
  double t1448;
  double t1508;
  double t1509;
  double t1512;
  double t1533;
  double t1534;
  double t1491;
  double t1501;
  double t1504;
  double t1505;
  double t1513;
  double t1514;
  double t1541;
  double t1543;
  double t1544;
  double t1545;
  double t1535;
  double t1536;
  double t1537;
  double t1538;
  double t1518;
  double t1519;
  double t1520;
  double t1522;
  double t1523;
  double t1524;
  double t1564;
  double t1563;
  double t1566;
  double t1559;
  double t1560;
  double t1592;
  double t1593;
  double t1604;
  double t1605;
  double t1606;
  double t1608;
  double t1609;
  double t1596;
  double t1597;
  double t1598;
  double t1599;
  double t1600;
  double t1549;
  double t1550;
  double t1551;
  double t1552;
  double t1553;
  double t1532;
  double t1539;
  double t1540;
  double t1546;
  double t1547;
  double t1627;
  double t1628;
  double t1629;
  double t1630;
  double t1631;
  double t1573;
  double t1574;
  double t1581;
  double t1583;
  double t1584;
  double t1616;
  double t1617;
  double t1618;
  double t1601;
  double t1610;
  double t1614;
  double t1661;
  double t1662;
  double t1663;
  double t1664;
  double t1665;
  double t1658;
  double t1659;
  double t1626;
  double t1632;
  double t1634;
  double t1635;
  double t1636;
  double t1684;
  double t1685;
  double t1686;
  double t1681;
  double t1682;
  double t1637;
  double t1641;
  double t1642;
  double t1643;
  double t1644;
  double t1645;
  double t1646;
  double t1647;
  double t1648;
  double t1649;
  double t1651;
  double t1653;
  double t1655;
  double t1656;
  double t1657;
  double t1667;
  double t1670;
  double t1672;
  double t1704;
  double t1705;
  double t1706;
  double t1700;
  double t1701;
  double t1702;
  double t1674;
  t1494 = Cos(var1[5]);
  t1472 = Cos(var1[6]);
  t1483 = Sin(var1[5]);
  t1495 = Sin(var1[6]);
  t1507 = Cos(var1[2]);
  t1448 = Sin(var1[2]);
  t1508 = t1494*t1472;
  t1509 = -1.*t1483*t1495;
  t1512 = t1508 + t1509;
  t1533 = -1.*t1472;
  t1534 = 1. + t1533;
  t1491 = -1.*t1472*t1483;
  t1501 = -1.*t1494*t1495;
  t1504 = t1491 + t1501;
  t1505 = -1.*t1448*t1504;
  t1513 = -1.*t1507*t1512;
  t1514 = t1505 + t1513;
  t1541 = -0.0265*t1534;
  t1543 = -0.025226*t1472;
  t1544 = -0.07700600000000002*t1495;
  t1545 = t1541 + t1543 + t1544;
  t1535 = -0.2375*t1534;
  t1536 = -0.314506*t1472;
  t1537 = -0.0012740000000000008*t1495;
  t1538 = t1535 + t1536 + t1537;
  t1518 = t1472*t1483;
  t1519 = t1494*t1495;
  t1520 = t1518 + t1519;
  t1522 = -1.*t1507*t1520;
  t1523 = -1.*t1448*t1512;
  t1524 = t1522 + t1523;
  t1564 = t1507*t1512;
  t1563 = -1.*t1448*t1520;
  t1566 = t1563 + t1564;
  t1559 = t1507*t1504;
  t1560 = t1559 + t1523;
  t1592 = -1.*t1494;
  t1593 = 1. + t1592;
  t1604 = -0.0695*t1593;
  t1605 = -0.0265*t1483;
  t1606 = -1.*t1483*t1545;
  t1608 = t1494*t1538;
  t1609 = t1604 + t1605 + t1606 + t1608;
  t1596 = -0.0265*t1593;
  t1597 = 0.0695*t1483;
  t1598 = t1494*t1545;
  t1599 = t1483*t1538;
  t1600 = t1596 + t1597 + t1598 + t1599;
  t1549 = -0.0265*t1472;
  t1550 = -1.*t1472*t1545;
  t1551 = 0.0695*t1495;
  t1552 = t1538*t1495;
  t1553 = t1549 + t1550 + t1551 + t1552;
  t1532 = 0.0695*t1472;
  t1539 = t1472*t1538;
  t1540 = 0.0265*t1495;
  t1546 = t1545*t1495;
  t1547 = t1532 + t1539 + t1540 + t1546;
  t1627 = -1.*t1494*t1472;
  t1628 = t1483*t1495;
  t1629 = t1627 + t1628;
  t1630 = t1507*t1629;
  t1631 = t1505 + t1630;
  t1573 = t1448*t1504;
  t1574 = t1573 + t1564;
  t1581 = t1507*t1520;
  t1583 = t1448*t1512;
  t1584 = t1581 + t1583;
  t1616 = -1.*t1609*t1504;
  t1617 = -1.*t1600*t1512;
  t1618 = t1616 + t1617;
  t1601 = t1600*t1520;
  t1610 = t1609*t1512;
  t1614 = t1601 + t1610;
  t1661 = -0.0265*t1494;
  t1662 = -0.0695*t1483;
  t1663 = -1.*t1494*t1545;
  t1664 = -1.*t1483*t1538;
  t1665 = t1661 + t1662 + t1663 + t1664;
  t1658 = 0.0695*t1494;
  t1659 = t1658 + t1605 + t1606 + t1608;
  t1626 = -0.0002543413600000002*t1560;
  t1632 = -0.015373477840000005*t1631;
  t1634 = t1626 + t1632;
  t1635 = 0.5*var2[6]*t1634;
  t1636 = 0.19964*t1553*t1560;
  t1684 = -0.0012740000000000008*t1472;
  t1685 = 0.07700600000000002*t1495;
  t1686 = t1684 + t1685;
  t1681 = -0.07700600000000002*t1472;
  t1682 = t1681 + t1537;
  t1637 = 0.19964*t1547*t1631;
  t1641 = 0.39928*t1566*t1560;
  t1642 = 0.39928*t1560*t1631;
  t1643 = t1641 + t1642;
  t1644 = 0.5*var2[1]*t1643;
  t1645 = 0.19964*t1574*t1566;
  t1646 = 0.19964*t1560*t1584;
  t1647 = 0.19964*t1574*t1631;
  t1648 = t1448*t1629;
  t1649 = t1559 + t1648;
  t1651 = 0.19964*t1560*t1649;
  t1653 = t1645 + t1646 + t1647 + t1651;
  t1655 = 0.5*var2[0]*t1653;
  t1656 = 0.19964*t1560*t1618;
  t1657 = t1609*t1504;
  t1667 = t1600*t1512;
  t1670 = 0.19964*t1614*t1631;
  t1672 = -1.*t1600*t1504;
  t1704 = -1.*t1483*t1682;
  t1705 = t1494*t1686;
  t1706 = t1704 + t1705;
  t1700 = t1494*t1682;
  t1701 = t1483*t1686;
  t1702 = t1700 + t1701;
  t1674 = -1.*t1609*t1629;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(0.5*(0.19964*Power(t1560,2) + 0.19964*Power(t1566,2) + 0.19964*t1514*t1574 + 0.19964*t1524*t1584)*var2[0] + 0.5*(0.39928*t1514*t1560 + 0.39928*t1524*t1566)*var2[1] + 0.5*(0.19964*t1514*t1614 + 0.19964*t1524*t1618)*var2[2] + 0.5*(0.19964*t1514*t1547 + 0.19964*t1524*t1553)*var2[5] + 0.5*(-0.015373477840000005*t1514 - 0.0002543413600000002*t1524)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(t1635 + t1644 + t1655 + 0.5*(t1656 + 0.19964*t1560*(t1657 + t1520*t1659 + t1512*t1665 + t1667) + t1670 + 0.19964*t1566*(-1.*t1512*t1659 - 1.*t1504*t1665 + t1672 + t1674))*var2[2] + 0.5*(t1636 + t1637)*var2[5]);
  p_output1[6]=var2[1]*(t1635 + t1644 + t1655 + 0.5*(t1656 + t1670 + 0.19964*t1566*(t1672 + t1674 - 1.*t1512*t1702 - 1.*t1504*t1706) + 0.19964*t1560*(t1657 + t1667 + t1520*t1702 + t1512*t1706))*var2[2] + 0.5*(t1636 + t1637 + 0.19964*t1560*(0.0265*t1472 - 0.0695*t1495 - 1.*t1495*t1538 + t1472*t1545 + t1495*t1682 + t1472*t1686) + 0.19964*t1566*(t1532 + t1539 + t1540 + t1546 - 1.*t1472*t1682 + t1495*t1686))*var2[5]);
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

#include "Ce3_vec_L5_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L5_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
