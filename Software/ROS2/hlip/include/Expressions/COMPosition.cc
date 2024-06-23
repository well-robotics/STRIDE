/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:36 GMT-05:00
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
static void output1(double *p_output1,const double *var1)
{
  double t1558;
  double t1611;
  double t1674;
  double t1684;
  double t1688;
  double t1691;
  double t1689;
  double t1693;
  double t1703;
  double t1708;
  double t1711;
  double t1716;
  double t1722;
  double t1725;
  double t1738;
  double t1739;
  double t1740;
  double t1726;
  double t1727;
  double t1729;
  double t1754;
  double t1755;
  double t1756;
  double t1760;
  double t1792;
  double t1793;
  double t1796;
  double t1798;
  double t1797;
  double t1801;
  double t1803;
  double t1804;
  double t1806;
  double t1807;
  double t1808;
  double t1813;
  double t1822;
  double t1825;
  double t1826;
  double t1828;
  double t1829;
  double t1831;
  double t1837;
  double t1840;
  double t1841;
  double t1844;
  double t1874;
  double t1875;
  double t1876;
  double t1877;
  double t1879;
  double t1757;
  double t1761;
  double t1762;
  double t1769;
  double t1770;
  double t1774;
  double t1776;
  double t1904;
  double t1905;
  double t1842;
  double t1845;
  double t1846;
  double t1907;
  double t1909;
  double t1910;
  double t1848;
  double t1850;
  double t1851;
  double t1853;
  t1558 = Cos(var1[2]);
  t1611 = Sin(var1[2]);
  t1674 = Cos(var1[3]);
  t1684 = -1.*t1674;
  t1688 = 1. + t1684;
  t1691 = Sin(var1[3]);
  t1689 = -0.0265*t1688;
  t1693 = -0.0695*t1691;
  t1703 = t1689 + t1693;
  t1708 = t1558*t1703;
  t1711 = -0.0695*t1688;
  t1716 = 0.0265*t1691;
  t1722 = t1711 + t1716;
  t1725 = t1611*t1722;
  t1738 = t1558*t1674;
  t1739 = t1611*t1691;
  t1740 = t1738 + t1739;
  t1726 = t1674*t1611;
  t1727 = -1.*t1558*t1691;
  t1729 = t1726 + t1727;
  t1754 = Cos(var1[4]);
  t1755 = -1.*t1754;
  t1756 = 1. + t1755;
  t1760 = Sin(var1[4]);
  t1792 = Cos(var1[5]);
  t1793 = -1.*t1792;
  t1796 = 1. + t1793;
  t1798 = Sin(var1[5]);
  t1797 = -0.0695*t1796;
  t1801 = -0.0265*t1798;
  t1803 = t1797 + t1801;
  t1804 = t1611*t1803;
  t1806 = -0.0265*t1796;
  t1807 = 0.0695*t1798;
  t1808 = t1806 + t1807;
  t1813 = t1558*t1808;
  t1822 = t1792*t1611;
  t1825 = t1558*t1798;
  t1826 = t1822 + t1825;
  t1828 = t1558*t1792;
  t1829 = -1.*t1611*t1798;
  t1831 = t1828 + t1829;
  t1837 = Cos(var1[6]);
  t1840 = -1.*t1837;
  t1841 = 1. + t1840;
  t1844 = Sin(var1[6]);
  t1874 = -1.*t1611*t1703;
  t1875 = t1558*t1722;
  t1876 = -1.*t1674*t1611;
  t1877 = t1558*t1691;
  t1879 = t1876 + t1877;
  t1757 = -0.0265*t1756;
  t1761 = -0.2375*t1760;
  t1762 = t1757 + t1761;
  t1769 = -0.2375*t1756;
  t1770 = 0.0265*t1760;
  t1774 = t1769 + t1770;
  t1776 = t1754*t1740;
  t1904 = t1558*t1803;
  t1905 = -1.*t1611*t1808;
  t1842 = -0.2375*t1841;
  t1845 = -0.0265*t1844;
  t1846 = t1842 + t1845;
  t1907 = -1.*t1792*t1611;
  t1909 = -1.*t1558*t1798;
  t1910 = t1907 + t1909;
  t1848 = -0.0265*t1841;
  t1850 = 0.2375*t1844;
  t1851 = t1848 + t1850;
  t1853 = t1837*t1831;
  p_output1[0]=0.28330859104971495*(1.7566*(-0.026461*t1558 - 0.046589*t1611 + var1[0]) + 0.69051*(t1708 + t1725 - 0.15232*t1729 - 0.025367*t1740 + var1[0]) + 0.19605*(t1708 + t1725 - 0.312342*(t1729*t1754 - 1.*t1740*t1760) + t1740*t1762 + t1729*t1774 - 0.025159*(t1729*t1760 + t1776) + var1[0]) + 0.69051*(t1804 + t1813 - 0.15232*t1826 - 0.025413*t1831 + var1[0]) + 0.19605*(t1804 + t1813 - 0.312334*(t1826*t1837 + t1831*t1844) + t1826*t1846 + t1831*t1851 - 0.025157*(-1.*t1826*t1844 + t1853) + var1[0]));
  p_output1[1]=0.0021663594534410656;
  p_output1[2]=0.28330859104971495*(1.7566*(-0.046589*t1558 + 0.026461*t1611 + var1[1]) + 0.69051*(-0.15232*t1740 + t1874 + t1875 - 0.025367*t1879 + var1[1]) + 0.19605*(t1740*t1774 + t1874 + t1875 + t1762*t1879 - 0.025159*(t1740*t1760 + t1754*t1879) - 0.312342*(t1776 - 1.*t1760*t1879) + var1[1]) + 0.69051*(-0.15232*t1831 + t1904 + t1905 - 0.025413*t1910 + var1[1]) + 0.19605*(t1831*t1846 + t1904 + t1905 + t1851*t1910 - 0.025157*(-1.*t1831*t1844 + t1837*t1910) - 0.312334*(t1853 + t1844*t1910) + var1[1]));
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

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
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

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 3, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "COMPosition.hh"

namespace SymFunction
{

void COMPosition_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
