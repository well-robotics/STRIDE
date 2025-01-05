/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:20 GMT-05:00
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
static void output1(double *p_output1,const double *var1)
{
  double t1615;
  double t1668;
  double t1731;
  double t1741;
  double t1745;
  double t1748;
  double t1746;
  double t1750;
  double t1760;
  double t1765;
  double t1768;
  double t1773;
  double t1779;
  double t1782;
  double t1795;
  double t1796;
  double t1797;
  double t1783;
  double t1784;
  double t1786;
  double t1811;
  double t1812;
  double t1813;
  double t1817;
  double t1849;
  double t1850;
  double t1853;
  double t1855;
  double t1854;
  double t1858;
  double t1860;
  double t1861;
  double t1863;
  double t1864;
  double t1865;
  double t1870;
  double t1879;
  double t1882;
  double t1883;
  double t1885;
  double t1886;
  double t1888;
  double t1894;
  double t1897;
  double t1898;
  double t1901;
  double t1931;
  double t1932;
  double t1933;
  double t1934;
  double t1936;
  double t1814;
  double t1818;
  double t1819;
  double t1826;
  double t1827;
  double t1831;
  double t1833;
  double t1961;
  double t1962;
  double t1899;
  double t1902;
  double t1903;
  double t1964;
  double t1966;
  double t1967;
  double t1905;
  double t1907;
  double t1908;
  double t1910;
  t1615 = Cos(var1[2]);
  t1668 = Sin(var1[2]);
  t1731 = Cos(var1[3]);
  t1741 = -1.*t1731;
  t1745 = 1. + t1741;
  t1748 = Sin(var1[3]);
  t1746 = -0.0265*t1745;
  t1750 = -0.0695*t1748;
  t1760 = t1746 + t1750;
  t1765 = t1615*t1760;
  t1768 = -0.0695*t1745;
  t1773 = 0.0265*t1748;
  t1779 = t1768 + t1773;
  t1782 = t1668*t1779;
  t1795 = t1615*t1731;
  t1796 = t1668*t1748;
  t1797 = t1795 + t1796;
  t1783 = t1731*t1668;
  t1784 = -1.*t1615*t1748;
  t1786 = t1783 + t1784;
  t1811 = Cos(var1[4]);
  t1812 = -1.*t1811;
  t1813 = 1. + t1812;
  t1817 = Sin(var1[4]);
  t1849 = Cos(var1[5]);
  t1850 = -1.*t1849;
  t1853 = 1. + t1850;
  t1855 = Sin(var1[5]);
  t1854 = -0.0695*t1853;
  t1858 = -0.0265*t1855;
  t1860 = t1854 + t1858;
  t1861 = t1668*t1860;
  t1863 = -0.0265*t1853;
  t1864 = 0.0695*t1855;
  t1865 = t1863 + t1864;
  t1870 = t1615*t1865;
  t1879 = t1849*t1668;
  t1882 = t1615*t1855;
  t1883 = t1879 + t1882;
  t1885 = t1615*t1849;
  t1886 = -1.*t1668*t1855;
  t1888 = t1885 + t1886;
  t1894 = Cos(var1[6]);
  t1897 = -1.*t1894;
  t1898 = 1. + t1897;
  t1901 = Sin(var1[6]);
  t1931 = -1.*t1668*t1760;
  t1932 = t1615*t1779;
  t1933 = -1.*t1731*t1668;
  t1934 = t1615*t1748;
  t1936 = t1933 + t1934;
  t1814 = -0.0265*t1813;
  t1818 = -0.2375*t1817;
  t1819 = t1814 + t1818;
  t1826 = -0.2375*t1813;
  t1827 = 0.0265*t1817;
  t1831 = t1826 + t1827;
  t1833 = t1811*t1797;
  t1961 = t1615*t1860;
  t1962 = -1.*t1668*t1865;
  t1899 = -0.2375*t1898;
  t1902 = -0.0265*t1901;
  t1903 = t1899 + t1902;
  t1964 = -1.*t1849*t1668;
  t1966 = -1.*t1615*t1855;
  t1967 = t1964 + t1966;
  t1905 = -0.0265*t1898;
  t1907 = 0.2375*t1901;
  t1908 = t1905 + t1907;
  t1910 = t1894*t1888;
  p_output1[0]=0.2996793431028799*(1.5566*(-0.026461*t1615 - 0.046589*t1668 + var1[0]) + 0.69051*(t1765 + t1782 - 0.15232*t1786 - 0.025367*t1797 + var1[0]) + 0.19964*(t1765 + t1782 - 0.314514*(t1786*t1811 - 1.*t1797*t1817) + t1797*t1819 + t1786*t1831 - 0.025229*(t1786*t1817 + t1833) + var1[0]) + 0.69051*(t1861 + t1870 - 0.15232*t1883 - 0.025413*t1888 + var1[0]) + 0.19964*(t1861 + t1870 - 0.314506*(t1883*t1894 + t1888*t1901) + t1883*t1903 + t1888*t1908 - 0.025226*(-1.*t1883*t1901 + t1910) + var1[0]));
  p_output1[1]=0.002039417585183853;
  p_output1[2]=0.2996793431028799*(1.5566*(-0.046589*t1615 + 0.026461*t1668 + var1[1]) + 0.69051*(-0.15232*t1797 + t1931 + t1932 - 0.025367*t1936 + var1[1]) + 0.19964*(t1797*t1831 + t1931 + t1932 + t1819*t1936 - 0.025229*(t1797*t1817 + t1811*t1936) - 0.314514*(t1833 - 1.*t1817*t1936) + var1[1]) + 0.69051*(-0.15232*t1888 + t1961 + t1962 - 0.025413*t1967 + var1[1]) + 0.19964*(t1888*t1903 + t1961 + t1962 + t1908*t1967 - 0.025226*(-1.*t1888*t1901 + t1894*t1967) - 0.314506*(t1910 + t1901*t1967) + var1[1]));
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
