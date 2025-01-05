/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:34:17 GMT-05:00
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
  double t1531;
  double t1557;
  double t1544;
  double t1549;
  double t1585;
  double t1581;
  double t1566;
  double t1569;
  double t1573;
  double t1574;
  double t1551;
  double t1554;
  double t1558;
  double t1560;
  double t1601;
  double t1578;
  double t1604;
  double t1605;
  double t1608;
  double t1619;
  double t1621;
  double t1623;
  double t1624;
  double t1625;
  double t1626;
  double t1628;
  double t1630;
  double t1536;
  double t1562;
  double t1564;
  double t1576;
  double t1577;
  double t1584;
  double t1590;
  double t1596;
  double t1639;
  double t1631;
  double t1640;
  double t1649;
  double t1614;
  double t1684;
  double t1685;
  double t1687;
  double t1707;
  double t1708;
  double t1709;
  double t1701;
  double t1703;
  double t1715;
  double t1716;
  double t1717;
  double t1718;
  double t1719;
  double t1722;
  double t1723;
  double t1725;
  double t1727;
  double t1729;
  double t1712;
  double t1713;
  double t1705;
  double t1710;
  double t1754;
  double t1755;
  double t1748;
  double t1749;
  double t1750;
  double t1746;
  double t1747;
  double t1751;
  double t1752;
  double t1753;
  double t1756;
  double t1757;
  double t1759;
  double t1760;
  double t1761;
  double t1597;
  double t1615;
  double t1683;
  double t1688;
  double t1689;
  double t1690;
  double t1693;
  double t1694;
  double t1695;
  double t1696;
  double t1697;
  double t1698;
  double t1711;
  double t1731;
  double t1737;
  double t1792;
  double t1793;
  double t1794;
  double t1788;
  double t1789;
  double t1790;
  double t1739;
  t1531 = Cos(var1[6]);
  t1557 = Sin(var1[6]);
  t1544 = -1.*t1531;
  t1549 = 1. + t1544;
  t1585 = Cos(var1[5]);
  t1581 = Sin(var1[5]);
  t1566 = -0.2375*t1549;
  t1569 = -0.314506*t1531;
  t1573 = -0.0012740000000000008*t1557;
  t1574 = t1566 + t1569 + t1573;
  t1551 = -0.0265*t1549;
  t1554 = -0.025226*t1531;
  t1558 = -0.07700600000000002*t1557;
  t1560 = t1551 + t1554 + t1558;
  t1601 = Cos(var1[2]);
  t1578 = Sin(var1[2]);
  t1604 = t1585*t1531;
  t1605 = -1.*t1581*t1557;
  t1608 = t1604 + t1605;
  t1619 = 0.0695*t1531;
  t1621 = t1531*t1574;
  t1623 = 0.0265*t1557;
  t1624 = t1560*t1557;
  t1625 = t1619 + t1621 + t1623 + t1624;
  t1626 = -1.*t1531*t1581;
  t1628 = -1.*t1585*t1557;
  t1630 = t1626 + t1628;
  t1536 = -0.0265*t1531;
  t1562 = -1.*t1531*t1560;
  t1564 = 0.0695*t1557;
  t1576 = t1574*t1557;
  t1577 = t1536 + t1562 + t1564 + t1576;
  t1584 = t1531*t1581;
  t1590 = t1585*t1557;
  t1596 = t1584 + t1590;
  t1639 = -1.*t1578*t1608;
  t1631 = t1601*t1630;
  t1640 = t1631 + t1639;
  t1649 = -1.*t1578*t1630;
  t1614 = t1601*t1608;
  t1684 = -1.*t1585*t1531;
  t1685 = t1581*t1557;
  t1687 = t1684 + t1685;
  t1707 = -0.0265*t1581;
  t1708 = -1.*t1581*t1560;
  t1709 = t1585*t1574;
  t1701 = -1.*t1585;
  t1703 = 1. + t1701;
  t1715 = -0.0265*t1585;
  t1716 = -0.0695*t1581;
  t1717 = -1.*t1585*t1560;
  t1718 = -1.*t1581*t1574;
  t1719 = t1715 + t1716 + t1717 + t1718;
  t1722 = -0.0265*t1703;
  t1723 = 0.0695*t1581;
  t1725 = t1585*t1560;
  t1727 = t1581*t1574;
  t1729 = t1722 + t1723 + t1725 + t1727;
  t1712 = 0.0695*t1585;
  t1713 = t1712 + t1707 + t1708 + t1709;
  t1705 = -0.0695*t1703;
  t1710 = t1705 + t1707 + t1708 + t1709;
  t1754 = -0.07700600000000002*t1531;
  t1755 = t1754 + t1573;
  t1748 = -0.0012740000000000008*t1531;
  t1749 = 0.07700600000000002*t1557;
  t1750 = t1748 + t1749;
  t1746 = 0.0265*t1531;
  t1747 = t1531*t1560;
  t1751 = t1531*t1750;
  t1752 = -0.0695*t1557;
  t1753 = -1.*t1574*t1557;
  t1756 = t1755*t1557;
  t1757 = t1746 + t1747 + t1751 + t1752 + t1753 + t1756;
  t1759 = -1.*t1531*t1755;
  t1760 = t1750*t1557;
  t1761 = t1619 + t1621 + t1759 + t1623 + t1624 + t1760;
  t1597 = -1.*t1578*t1596;
  t1615 = t1597 + t1614;
  t1683 = 0.19964*t1577*t1640;
  t1688 = t1601*t1687;
  t1689 = t1649 + t1688;
  t1690 = 0.19964*t1625*t1689;
  t1693 = t1578*t1630;
  t1694 = t1693 + t1614;
  t1695 = 0.19964*t1577*t1694;
  t1696 = t1578*t1687;
  t1697 = t1631 + t1696;
  t1698 = 0.19964*t1625*t1697;
  t1711 = t1710*t1630;
  t1731 = t1729*t1608;
  t1737 = -1.*t1729*t1630;
  t1792 = -1.*t1581*t1755;
  t1793 = t1585*t1750;
  t1794 = t1792 + t1793;
  t1788 = t1585*t1755;
  t1789 = t1581*t1750;
  t1790 = t1788 + t1789;
  t1739 = -1.*t1710*t1687;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.19964*t1577*t1615 + 0.19964*t1625*t1640)*var2[0] + 0.5*(0.19964*t1577*(-1.*t1596*t1601 + t1639) + 0.19964*t1625*(-1.*t1601*t1608 + t1649))*var2[1])*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(0.5*(t1695 + t1698)*var2[0] + 0.5*(t1683 + t1690)*var2[1] + 0.5*(0.19964*t1625*(t1711 + t1596*t1713 + t1608*t1719 + t1731) + 0.19964*t1577*(-1.*t1608*t1713 - 1.*t1630*t1719 + t1737 + t1739))*var2[2])*var2[5];
  p_output1[6]=var2[5]*(0.5*(t1695 + t1698 + 0.19964*t1694*t1757 + 0.19964*(t1596*t1601 + t1578*t1608)*t1761)*var2[0] + 0.5*(t1683 + t1690 + 0.19964*t1640*t1757 + 0.19964*t1615*t1761)*var2[1] + 0.5*(0.19964*(t1608*t1710 + t1596*t1729)*t1757 + 0.19964*(-1.*t1630*t1710 - 1.*t1608*t1729)*t1761 + 0.19964*t1625*(t1711 + t1731 + t1596*t1790 + t1608*t1794) + 0.19964*t1577*(t1737 + t1739 - 1.*t1608*t1790 - 1.*t1630*t1794))*var2[2] + 0.5*(0.39928*t1625*t1757 + 0.39928*t1577*t1761)*var2[5] + 0.5*(-0.015373477840000005*t1757 - 0.0002543413600000002*t1761)*var2[6]);
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

#include "Ce3_vec_L5_J6_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L5_J6_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
