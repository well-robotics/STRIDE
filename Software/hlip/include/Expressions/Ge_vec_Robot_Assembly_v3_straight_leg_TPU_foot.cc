/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:27 GMT-05:00
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
  double t1501;
  double t1519;
  double t1567;
  double t1574;
  double t1583;
  double t1590;
  double t1584;
  double t1594;
  double t1604;
  double t1610;
  double t1611;
  double t1612;
  double t1613;
  double t1620;
  double t1632;
  double t1635;
  double t1636;
  double t1621;
  double t1622;
  double t1624;
  double t1648;
  double t1650;
  double t1652;
  double t1664;
  double t1704;
  double t1706;
  double t1707;
  double t1709;
  double t1708;
  double t1711;
  double t1713;
  double t1714;
  double t1715;
  double t1716;
  double t1717;
  double t1718;
  double t1719;
  double t1720;
  double t1721;
  double t1723;
  double t1724;
  double t1725;
  double t1730;
  double t1731;
  double t1732;
  double t1736;
  double t1759;
  double t1760;
  double t1761;
  double t1762;
  double t1763;
  double t1764;
  double t1765;
  double t1770;
  double t1771;
  double t1772;
  double t1654;
  double t1665;
  double t1670;
  double t1766;
  double t1767;
  double t1768;
  double t1675;
  double t1681;
  double t1683;
  double t1778;
  double t1783;
  double t1722;
  double t1726;
  double t1804;
  double t1805;
  double t1806;
  double t1807;
  double t1808;
  double t1809;
  double t1810;
  double t1734;
  double t1738;
  double t1739;
  double t1740;
  double t1741;
  double t1742;
  double t1743;
  double t1744;
  double t1745;
  double t1746;
  double t1748;
  double t1750;
  double t1752;
  double t1753;
  double t1754;
  double t1755;
  double t1816;
  double t1817;
  double t1818;
  t1501 = Cos(var1[2]);
  t1519 = Sin(var1[2]);
  t1567 = Cos(var1[3]);
  t1574 = -1.*t1567;
  t1583 = 1. + t1574;
  t1590 = Sin(var1[3]);
  t1584 = -0.0265*t1583;
  t1594 = -0.0695*t1590;
  t1604 = t1584 + t1594;
  t1610 = -1.*t1501*t1604;
  t1611 = -0.0695*t1583;
  t1612 = 0.0265*t1590;
  t1613 = t1611 + t1612;
  t1620 = -1.*t1519*t1613;
  t1632 = -1.*t1501*t1567;
  t1635 = -1.*t1519*t1590;
  t1636 = t1632 + t1635;
  t1621 = -1.*t1567*t1519;
  t1622 = t1501*t1590;
  t1624 = t1621 + t1622;
  t1648 = Cos(var1[4]);
  t1650 = -1.*t1648;
  t1652 = 1. + t1650;
  t1664 = Sin(var1[4]);
  t1704 = Cos(var1[5]);
  t1706 = -1.*t1704;
  t1707 = 1. + t1706;
  t1709 = Sin(var1[5]);
  t1708 = -0.0695*t1707;
  t1711 = -0.0265*t1709;
  t1713 = t1708 + t1711;
  t1714 = -1.*t1519*t1713;
  t1715 = -0.0265*t1707;
  t1716 = 0.0695*t1709;
  t1717 = t1715 + t1716;
  t1718 = -1.*t1501*t1717;
  t1719 = -1.*t1704*t1519;
  t1720 = -1.*t1501*t1709;
  t1721 = t1719 + t1720;
  t1723 = -1.*t1501*t1704;
  t1724 = t1519*t1709;
  t1725 = t1723 + t1724;
  t1730 = Cos(var1[6]);
  t1731 = -1.*t1730;
  t1732 = 1. + t1731;
  t1736 = Sin(var1[6]);
  t1759 = 0.0265*t1567;
  t1760 = t1759 + t1594;
  t1761 = t1501*t1760;
  t1762 = -0.0695*t1567;
  t1763 = -0.0265*t1590;
  t1764 = t1762 + t1763;
  t1765 = -1.*t1519*t1764;
  t1770 = t1501*t1567;
  t1771 = t1519*t1590;
  t1772 = t1770 + t1771;
  t1654 = -0.0265*t1652;
  t1665 = -0.2375*t1664;
  t1670 = t1654 + t1665;
  t1766 = t1567*t1519;
  t1767 = -1.*t1501*t1590;
  t1768 = t1766 + t1767;
  t1675 = -0.2375*t1652;
  t1681 = 0.0265*t1664;
  t1683 = t1675 + t1681;
  t1778 = t1648*t1772;
  t1783 = -1.*t1772*t1664;
  t1722 = -0.15232*t1721;
  t1726 = -0.025413*t1725;
  t1804 = -0.0265*t1704;
  t1805 = -0.0695*t1709;
  t1806 = t1804 + t1805;
  t1807 = t1501*t1806;
  t1808 = 0.0695*t1704;
  t1809 = t1808 + t1711;
  t1810 = -1.*t1519*t1809;
  t1734 = -0.2375*t1732;
  t1738 = -0.0265*t1736;
  t1739 = t1734 + t1738;
  t1740 = t1721*t1739;
  t1741 = -0.0265*t1732;
  t1742 = 0.2375*t1736;
  t1743 = t1741 + t1742;
  t1744 = t1725*t1743;
  t1745 = t1730*t1725;
  t1746 = -1.*t1721*t1736;
  t1748 = t1745 + t1746;
  t1750 = -0.025157*t1748;
  t1752 = t1730*t1721;
  t1753 = t1725*t1736;
  t1754 = t1752 + t1753;
  t1755 = -0.312334*t1754;
  t1816 = t1501*t1704;
  t1817 = -1.*t1519*t1709;
  t1818 = t1816 + t1817;
  p_output1[0]=0;
  p_output1[1]=-34.626553200000004;
  p_output1[2]=-17.232246*(0.026461*t1501 + 0.046589*t1519) - 6.7739031*(t1610 + t1620 - 0.15232*t1624 - 0.025367*t1636) - 1.9232505000000002*(t1610 + t1620 - 0.025159*(t1636*t1648 + t1624*t1664) - 0.312342*(t1624*t1648 - 1.*t1636*t1664) + t1636*t1670 + t1624*t1683) - 6.7739031*(t1714 + t1718 + t1722 + t1726) - 1.9232505000000002*(t1714 + t1718 + t1740 + t1744 + t1750 + t1755);
  p_output1[3]=-6.7739031*(t1761 + t1765 - 0.15232*t1768 - 0.025367*t1772) - 1.9232505000000002*(t1761 + t1765 + t1683*t1768 + t1670*t1772 - 0.025159*(t1664*t1768 + t1778) - 0.312342*(t1648*t1768 + t1783));
  p_output1[4]=-1.9232505000000002*(t1624*(-0.2375*t1648 - 0.0265*t1664) + (0.0265*t1648 + t1665)*t1772 - 0.025159*(-1.*t1624*t1664 + t1778) - 0.312342*(-1.*t1624*t1648 + t1783));
  p_output1[5]=-6.7739031*(t1722 + t1726 + t1807 + t1810) - 1.9232505000000002*(t1740 + t1744 + t1750 + t1755 + t1807 + t1810);
  p_output1[6]=-1.9232505000000002*(t1721*(0.2375*t1730 + t1738) + (-0.0265*t1730 - 0.2375*t1736)*t1818 - 0.025157*(t1746 - 1.*t1730*t1818) - 0.312334*(t1752 - 1.*t1736*t1818));
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

#include "Ge_vec_Robot_Assembly_v3_straight_leg_TPU_foot.hh"

namespace SymFunction
{

void Ge_vec_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
