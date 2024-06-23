/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:36 GMT-05:00
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
  double t1500;
  double t1484;
  double t1492;
  double t1501;
  double t1562;
  double t1499;
  double t1502;
  double t1524;
  double t1448;
  double t1573;
  double t1574;
  double t1575;
  double t1582;
  double t1583;
  double t1590;
  double t1591;
  double t1592;
  double t1593;
  double t1594;
  double t1595;
  double t1601;
  double t1561;
  double t1563;
  double t1564;
  double t1569;
  double t1570;
  double t1571;
  double t1605;
  double t1606;
  double t1607;
  double t1608;
  double t1609;
  double t1610;
  double t1625;
  double t1626;
  double t1635;
  double t1636;
  double t1637;
  double t1639;
  double t1640;
  double t1641;
  double t1645;
  double t1646;
  double t1647;
  double t1651;
  double t1652;
  double t1628;
  double t1629;
  double t1630;
  double t1602;
  double t1603;
  double t1604;
  double t1622;
  double t1623;
  double t1584;
  double t1585;
  double t1589;
  double t1597;
  double t1598;
  double t1599;
  double t1612;
  double t1613;
  double t1614;
  double t1624;
  double t1627;
  double t1631;
  double t1632;
  double t1633;
  double t1638;
  double t1642;
  double t1643;
  double t1648;
  double t1649;
  double t1650;
  double t1653;
  double t1654;
  double t1656;
  double t1657;
  double t1658;
  double t1660;
  double t1661;
  double t1662;
  double t1663;
  double t1664;
  double t1682;
  double t1683;
  double t1684;
  double t1685;
  double t1686;
  double t1644;
  double t1655;
  double t1659;
  double t1665;
  double t1666;
  double t1671;
  double t1672;
  double t1673;
  double t1674;
  double t1675;
  double t1596;
  double t1600;
  double t1611;
  double t1615;
  double t1616;
  double t1691;
  double t1692;
  double t1693;
  double t1694;
  double t1695;
  t1500 = Cos(var1[5]);
  t1484 = Cos(var1[6]);
  t1492 = Sin(var1[5]);
  t1501 = Sin(var1[6]);
  t1562 = Sin(var1[2]);
  t1499 = -1.*t1484*t1492;
  t1502 = -1.*t1500*t1501;
  t1524 = t1499 + t1502;
  t1448 = Cos(var1[2]);
  t1573 = -1.*t1484;
  t1574 = 1. + t1573;
  t1575 = 0.5*t1574;
  t1582 = 0.671885*t1484;
  t1583 = t1575 + t1582;
  t1590 = t1562*t1524;
  t1591 = t1500*t1484;
  t1592 = -1.*t1492*t1501;
  t1593 = t1591 + t1592;
  t1594 = t1448*t1593;
  t1595 = t1590 + t1594;
  t1601 = t1583*t1484;
  t1561 = t1448*t1524;
  t1563 = -1.*t1500*t1484;
  t1564 = t1492*t1501;
  t1569 = t1563 + t1564;
  t1570 = t1562*t1569;
  t1571 = t1561 + t1570;
  t1605 = t1484*t1492;
  t1606 = t1500*t1501;
  t1607 = t1605 + t1606;
  t1608 = t1448*t1607;
  t1609 = t1562*t1593;
  t1610 = t1608 + t1609;
  t1625 = -1.*t1562*t1593;
  t1626 = t1561 + t1625;
  t1635 = t1583*t1492;
  t1636 = 0.171885*t1500*t1501;
  t1637 = t1635 + t1636;
  t1639 = t1500*t1583;
  t1640 = -0.171885*t1492*t1501;
  t1641 = t1639 + t1640;
  t1645 = -0.171885*t1484*t1492;
  t1646 = -0.171885*t1500*t1501;
  t1647 = t1645 + t1646;
  t1651 = 0.171885*t1500*t1484;
  t1652 = t1651 + t1640;
  t1628 = -1.*t1562*t1524;
  t1629 = t1448*t1569;
  t1630 = t1628 + t1629;
  t1602 = Power(t1484,2);
  t1603 = -0.171885*t1602;
  t1604 = t1601 + t1603;
  t1622 = -1.*t1562*t1607;
  t1623 = t1622 + t1594;
  t1584 = t1583*t1501;
  t1585 = -0.171885*t1484*t1501;
  t1589 = t1584 + t1585;
  t1597 = -1.*t1583*t1501;
  t1598 = 0.171885*t1484*t1501;
  t1599 = t1597 + t1598;
  t1612 = Power(t1501,2);
  t1613 = 0.171885*t1612;
  t1614 = t1601 + t1613;
  t1624 = 0.85216*t1595*t1623;
  t1627 = 0.85216*t1626*t1610;
  t1631 = 0.85216*t1595*t1630;
  t1632 = 0.85216*t1626*t1571;
  t1633 = t1624 + t1627 + t1631 + t1632;
  t1638 = -1.*t1637*t1593;
  t1642 = -1.*t1524*t1641;
  t1643 = t1638 + t1642;
  t1648 = t1647*t1593;
  t1649 = t1637*t1593;
  t1650 = t1524*t1641;
  t1653 = t1607*t1652;
  t1654 = t1648 + t1649 + t1650 + t1653;
  t1656 = t1637*t1607;
  t1657 = t1593*t1641;
  t1658 = t1656 + t1657;
  t1660 = -1.*t1524*t1647;
  t1661 = -1.*t1524*t1637;
  t1662 = -1.*t1593*t1652;
  t1663 = -1.*t1641*t1569;
  t1664 = t1660 + t1661 + t1662 + t1663;
  t1682 = 0.85216*t1626*t1643;
  t1683 = 0.85216*t1626*t1654;
  t1684 = 0.85216*t1658*t1630;
  t1685 = 0.85216*t1623*t1664;
  t1686 = t1682 + t1683 + t1684 + t1685;
  t1644 = 0.85216*t1595*t1643;
  t1655 = 0.85216*t1595*t1654;
  t1659 = 0.85216*t1658*t1571;
  t1665 = 0.85216*t1610*t1664;
  t1666 = t1644 + t1655 + t1659 + t1665;
  t1671 = 0.85216*t1604*t1623;
  t1672 = 0.85216*t1589*t1626;
  t1673 = 0.85216*t1599*t1626;
  t1674 = 0.85216*t1614*t1630;
  t1675 = t1671 + t1672 + t1673 + t1674;
  t1596 = 0.85216*t1589*t1595;
  t1600 = 0.85216*t1599*t1595;
  t1611 = 0.85216*t1604*t1610;
  t1615 = 0.85216*t1614*t1571;
  t1616 = t1596 + t1600 + t1611 + t1615;
  t1691 = 0.85216*t1604*t1643;
  t1692 = 0.85216*t1599*t1658;
  t1693 = 0.85216*t1614*t1654;
  t1694 = 0.85216*t1589*t1664;
  t1695 = t1691 + t1692 + t1693 + t1694;
  p_output1[0]=var2[6]*(-0.5*(1.70432*t1571*t1595 + 1.70432*t1595*t1610)*var2[0] - 0.5*t1633*var2[1] - 0.5*t1666*var2[2] - 0.5*t1616*var2[5] - 0.0732367608*t1571*var2[6]);
  p_output1[1]=var2[6]*(-0.5*t1633*var2[0] - 0.5*(1.70432*t1623*t1626 + 1.70432*t1626*t1630)*var2[1] - 0.5*t1686*var2[2] - 0.5*t1675*var2[5] - 0.0732367608*t1630*var2[6]);
  p_output1[2]=var2[6]*(-0.5*t1666*var2[0] - 0.5*t1686*var2[1] - 0.5*(1.70432*t1654*t1658 + 1.70432*t1643*t1664)*var2[2] - 0.5*t1695*var2[5] - 0.0732367608*t1654*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[6]*(-0.5*t1616*var2[0] - 0.5*t1675*var2[1] - 0.5*t1695*var2[2] - 0.5*(1.70432*t1589*t1604 + 1.70432*t1599*t1614)*var2[5] - 0.0732367608*t1599*var2[6]);
  p_output1[6]=(-0.0732367608*t1571*var2[0] - 0.0732367608*t1630*var2[1] - 0.0732367608*t1654*var2[2] - 0.0732367608*t1599*var2[5])*var2[6];
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

#include "Ce2_vec7_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
