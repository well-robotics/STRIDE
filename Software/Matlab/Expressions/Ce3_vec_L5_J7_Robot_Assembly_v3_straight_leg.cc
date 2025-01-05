/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:34:21 GMT-05:00
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
  double t1564;
  double t1551;
  double t1554;
  double t1566;
  double t1577;
  double t1536;
  double t1584;
  double t1590;
  double t1604;
  double t1619;
  double t1621;
  double t1623;
  double t1558;
  double t1569;
  double t1573;
  double t1625;
  double t1624;
  double t1631;
  double t1648;
  double t1614;
  double t1688;
  double t1689;
  double t1691;
  double t1723;
  double t1725;
  double t1722;
  double t1727;
  double t1731;
  double t1732;
  double t1733;
  double t1736;
  double t1738;
  double t1740;
  double t1741;
  double t1742;
  double t1743;
  double t1715;
  double t1716;
  double t1749;
  double t1750;
  double t1751;
  double t1752;
  double t1754;
  double t1757;
  double t1758;
  double t1760;
  double t1761;
  double t1762;
  double t1746;
  double t1747;
  double t1721;
  double t1744;
  double t1784;
  double t1785;
  double t1778;
  double t1779;
  double t1780;
  double t1685;
  double t1692;
  double t1693;
  double t1694;
  double t1696;
  double t1697;
  double t1699;
  double t1700;
  double t1705;
  double t1707;
  double t1709;
  double t1711;
  double t1712;
  double t1714;
  double t1745;
  double t1763;
  double t1767;
  double t1806;
  double t1807;
  double t1808;
  double t1802;
  double t1803;
  double t1804;
  double t1769;
  t1564 = Cos(var1[5]);
  t1551 = Cos(var1[6]);
  t1554 = Sin(var1[5]);
  t1566 = Sin(var1[6]);
  t1577 = Cos(var1[2]);
  t1536 = Sin(var1[2]);
  t1584 = t1564*t1551;
  t1590 = -1.*t1554*t1566;
  t1604 = t1584 + t1590;
  t1619 = -1.*t1551*t1554;
  t1621 = -1.*t1564*t1566;
  t1623 = t1619 + t1621;
  t1558 = t1551*t1554;
  t1569 = t1564*t1566;
  t1573 = t1558 + t1569;
  t1625 = -1.*t1536*t1604;
  t1624 = t1577*t1623;
  t1631 = t1624 + t1625;
  t1648 = -1.*t1536*t1623;
  t1614 = t1577*t1604;
  t1688 = -1.*t1564*t1551;
  t1689 = t1554*t1566;
  t1691 = t1688 + t1689;
  t1723 = -1.*t1551;
  t1725 = 1. + t1723;
  t1722 = -0.0265*t1554;
  t1727 = -0.0265*t1725;
  t1731 = -0.025226*t1551;
  t1732 = -0.07700600000000002*t1566;
  t1733 = t1727 + t1731 + t1732;
  t1736 = -1.*t1554*t1733;
  t1738 = -0.2375*t1725;
  t1740 = -0.314506*t1551;
  t1741 = -0.0012740000000000008*t1566;
  t1742 = t1738 + t1740 + t1741;
  t1743 = t1564*t1742;
  t1715 = -1.*t1564;
  t1716 = 1. + t1715;
  t1749 = -0.0265*t1564;
  t1750 = -0.0695*t1554;
  t1751 = -1.*t1564*t1733;
  t1752 = -1.*t1554*t1742;
  t1754 = t1749 + t1750 + t1751 + t1752;
  t1757 = -0.0265*t1716;
  t1758 = 0.0695*t1554;
  t1760 = t1564*t1733;
  t1761 = t1554*t1742;
  t1762 = t1757 + t1758 + t1760 + t1761;
  t1746 = 0.0695*t1564;
  t1747 = t1746 + t1722 + t1736 + t1743;
  t1721 = -0.0695*t1716;
  t1744 = t1721 + t1722 + t1736 + t1743;
  t1784 = -0.07700600000000002*t1551;
  t1785 = t1784 + t1741;
  t1778 = -0.0012740000000000008*t1551;
  t1779 = 0.07700600000000002*t1566;
  t1780 = t1778 + t1779;
  t1685 = -0.0002543413600000002*t1631;
  t1692 = t1577*t1691;
  t1693 = t1648 + t1692;
  t1694 = -0.015373477840000005*t1693;
  t1696 = t1685 + t1694;
  t1697 = 0.5*var2[1]*t1696;
  t1699 = t1536*t1623;
  t1700 = t1699 + t1614;
  t1705 = -0.0002543413600000002*t1700;
  t1707 = t1536*t1691;
  t1709 = t1624 + t1707;
  t1711 = -0.015373477840000005*t1709;
  t1712 = t1705 + t1711;
  t1714 = 0.5*var2[0]*t1712;
  t1745 = t1744*t1623;
  t1763 = t1762*t1604;
  t1767 = -1.*t1762*t1623;
  t1806 = -1.*t1554*t1785;
  t1807 = t1564*t1780;
  t1808 = t1806 + t1807;
  t1802 = t1564*t1785;
  t1803 = t1554*t1780;
  t1804 = t1802 + t1803;
  t1769 = -1.*t1744*t1691;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(-0.0002543413600000002*(-1.*t1536*t1573 + t1614) - 0.015373477840000005*t1631)*var2[0] + 0.5*(-0.0002543413600000002*(-1.*t1573*t1577 + t1625) - 0.015373477840000005*(-1.*t1577*t1604 + t1648))*var2[1])*var2[6];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(t1697 + t1714 + 0.5*(-0.015373477840000005*(t1745 + t1573*t1747 + t1604*t1754 + t1763) - 0.0002543413600000002*(-1.*t1604*t1747 - 1.*t1623*t1754 + t1767 + t1769))*var2[2])*var2[6];
  p_output1[6]=(t1697 + t1714 + 0.5*(-0.015373477840000005*(t1745 + t1763 + t1573*t1804 + t1604*t1808) - 0.0002543413600000002*(t1767 + t1769 - 1.*t1604*t1804 - 1.*t1623*t1808))*var2[2] + 0.5*(-0.0002543413600000002*(0.0695*t1551 + 0.0265*t1566 + t1566*t1733 + t1551*t1742 + t1566*t1780 - 1.*t1551*t1785) - 0.015373477840000005*(0.0265*t1551 - 0.0695*t1566 + t1551*t1733 - 1.*t1566*t1742 + t1551*t1780 + t1566*t1785))*var2[5])*var2[6];
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

#include "Ce3_vec_L5_J7_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L5_J7_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
