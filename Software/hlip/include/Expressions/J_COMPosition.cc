/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:37 GMT-05:00
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
  double t1642;
  double t1584;
  double t1708;
  double t1711;
  double t1716;
  double t1725;
  double t1722;
  double t1726;
  double t1734;
  double t1738;
  double t1739;
  double t1741;
  double t1744;
  double t1745;
  double t1757;
  double t1761;
  double t1762;
  double t1769;
  double t1770;
  double t1774;
  double t1784;
  double t1785;
  double t1786;
  double t1788;
  double t1845;
  double t1846;
  double t1847;
  double t1850;
  double t1848;
  double t1851;
  double t1852;
  double t1855;
  double t1856;
  double t1857;
  double t1858;
  double t1860;
  double t1865;
  double t1866;
  double t1867;
  double t1861;
  double t1862;
  double t1863;
  double t1875;
  double t1877;
  double t1880;
  double t1883;
  double t1930;
  double t1931;
  double t1933;
  double t1934;
  double t1936;
  double t1787;
  double t1789;
  double t1791;
  double t1801;
  double t1803;
  double t1804;
  double t1827;
  double t1957;
  double t1962;
  double t1882;
  double t1884;
  double t1885;
  double t1969;
  double t1970;
  double t1976;
  double t1890;
  double t1893;
  double t1894;
  double t1902;
  double t1999;
  double t2000;
  double t2001;
  double t2002;
  double t2003;
  double t2004;
  double t2006;
  double t2007;
  double t2008;
  double t2012;
  double t1945;
  double t2042;
  double t2044;
  double t1807;
  double t2034;
  double t2030;
  double t2063;
  double t2065;
  double t2074;
  double t2075;
  double t2077;
  double t2078;
  double t2079;
  double t1813;
  double t1822;
  double t1864;
  double t1869;
  double t2107;
  double t2108;
  double t2109;
  double t2113;
  double t2115;
  double t2116;
  double t2118;
  double t1888;
  double t1895;
  double t1896;
  double t1897;
  double t1898;
  double t1901;
  double t1904;
  double t1911;
  double t1916;
  double t1966;
  double t1977;
  double t2129;
  double t2130;
  double t1980;
  double t1981;
  double t1982;
  double t1983;
  double t1984;
  double t1985;
  double t1987;
  double t1988;
  double t1989;
  double t2138;
  double t2139;
  double t2140;
  double t2142;
  double t2143;
  double t2144;
  double t2146;
  double t2147;
  t1642 = Sin(var1[2]);
  t1584 = Cos(var1[2]);
  t1708 = Cos(var1[3]);
  t1711 = -1.*t1708;
  t1716 = 1. + t1711;
  t1725 = Sin(var1[3]);
  t1722 = -0.0265*t1716;
  t1726 = -0.0695*t1725;
  t1734 = t1722 + t1726;
  t1738 = -1.*t1642*t1734;
  t1739 = -0.0695*t1716;
  t1741 = 0.0265*t1725;
  t1744 = t1739 + t1741;
  t1745 = t1584*t1744;
  t1757 = -1.*t1708*t1642;
  t1761 = t1584*t1725;
  t1762 = t1757 + t1761;
  t1769 = t1584*t1708;
  t1770 = t1642*t1725;
  t1774 = t1769 + t1770;
  t1784 = Cos(var1[4]);
  t1785 = -1.*t1784;
  t1786 = 1. + t1785;
  t1788 = Sin(var1[4]);
  t1845 = Cos(var1[5]);
  t1846 = -1.*t1845;
  t1847 = 1. + t1846;
  t1850 = Sin(var1[5]);
  t1848 = -0.0695*t1847;
  t1851 = -0.0265*t1850;
  t1852 = t1848 + t1851;
  t1855 = t1584*t1852;
  t1856 = -0.0265*t1847;
  t1857 = 0.0695*t1850;
  t1858 = t1856 + t1857;
  t1860 = -1.*t1642*t1858;
  t1865 = t1584*t1845;
  t1866 = -1.*t1642*t1850;
  t1867 = t1865 + t1866;
  t1861 = -1.*t1845*t1642;
  t1862 = -1.*t1584*t1850;
  t1863 = t1861 + t1862;
  t1875 = Cos(var1[6]);
  t1877 = -1.*t1875;
  t1880 = 1. + t1877;
  t1883 = Sin(var1[6]);
  t1930 = -1.*t1584*t1734;
  t1931 = -1.*t1642*t1744;
  t1933 = -1.*t1584*t1708;
  t1934 = -1.*t1642*t1725;
  t1936 = t1933 + t1934;
  t1787 = -0.0265*t1786;
  t1789 = -0.2375*t1788;
  t1791 = t1787 + t1789;
  t1801 = -0.2375*t1786;
  t1803 = 0.0265*t1788;
  t1804 = t1801 + t1803;
  t1827 = t1784*t1762;
  t1957 = -1.*t1642*t1852;
  t1962 = -1.*t1584*t1858;
  t1882 = -0.2375*t1880;
  t1884 = -0.0265*t1883;
  t1885 = t1882 + t1884;
  t1969 = -1.*t1584*t1845;
  t1970 = t1642*t1850;
  t1976 = t1969 + t1970;
  t1890 = -0.0265*t1880;
  t1893 = 0.2375*t1883;
  t1894 = t1890 + t1893;
  t1902 = t1875*t1863;
  t1999 = 0.0265*t1708;
  t2000 = t1999 + t1726;
  t2001 = t1642*t2000;
  t2002 = -0.0695*t1708;
  t2003 = -0.0265*t1725;
  t2004 = t2002 + t2003;
  t2006 = t1584*t2004;
  t2007 = t1708*t1642;
  t2008 = -1.*t1584*t1725;
  t2012 = t2007 + t2008;
  t1945 = t1784*t1936;
  t2042 = t1584*t2000;
  t2044 = -1.*t1642*t2004;
  t1807 = t1784*t1774;
  t2034 = t1784*t2012;
  t2030 = -1.*t2012*t1788;
  t2063 = -1.*t1774*t1788;
  t2065 = t2034 + t2063;
  t2074 = 0.0265*t1784;
  t2075 = t2074 + t1789;
  t2077 = -0.2375*t1784;
  t2078 = -0.0265*t1788;
  t2079 = t2077 + t2078;
  t1813 = -1.*t1762*t1788;
  t1822 = t1807 + t1813;
  t1864 = -0.025413*t1863;
  t1869 = -0.15232*t1867;
  t2107 = -0.0265*t1845;
  t2108 = -0.0695*t1850;
  t2109 = t2107 + t2108;
  t2113 = t1642*t2109;
  t2115 = 0.0695*t1845;
  t2116 = t2115 + t1851;
  t2118 = t1584*t2116;
  t1888 = t1867*t1885;
  t1895 = t1863*t1894;
  t1896 = t1875*t1867;
  t1897 = t1863*t1883;
  t1898 = t1896 + t1897;
  t1901 = -0.312334*t1898;
  t1904 = -1.*t1867*t1883;
  t1911 = t1902 + t1904;
  t1916 = -0.025157*t1911;
  t1966 = -0.15232*t1863;
  t1977 = -0.025413*t1976;
  t2129 = t1584*t2109;
  t2130 = -1.*t1642*t2116;
  t1980 = t1863*t1885;
  t1981 = t1976*t1894;
  t1982 = t1875*t1976;
  t1983 = -1.*t1863*t1883;
  t1984 = t1982 + t1983;
  t1985 = -0.025157*t1984;
  t1987 = t1976*t1883;
  t1988 = t1902 + t1987;
  t1989 = -0.312334*t1988;
  t2138 = t1845*t1642;
  t2139 = t1584*t1850;
  t2140 = t2138 + t2139;
  t2142 = -0.0265*t1875;
  t2143 = -0.2375*t1883;
  t2144 = t2142 + t2143;
  t2146 = 0.2375*t1875;
  t2147 = t2146 + t1884;
  p_output1[0]=0.9999999999999998;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=0.9999999999999998;
  p_output1[6]=0.28330859104971495*(1.7566*(-0.046589*t1584 + 0.026461*t1642) + 0.69051*(t1738 + t1745 - 0.025367*t1762 - 0.15232*t1774) + 0.19605*(t1738 + t1745 + t1762*t1791 + t1774*t1804 - 0.312342*t1822 - 0.025159*(t1774*t1788 + t1827)) + 0.69051*(t1855 + t1860 + t1864 + t1869) + 0.19605*(t1855 + t1860 + t1888 + t1895 + t1901 + t1916));
  p_output1[7]=0;
  p_output1[8]=0.28330859104971495*(1.7566*(0.026461*t1584 + 0.046589*t1642) + 0.69051*(-0.15232*t1762 + t1930 + t1931 - 0.025367*t1936) + 0.19605*(t1762*t1804 + t1930 + t1931 + t1791*t1936 - 0.312342*(t1827 - 1.*t1788*t1936) - 0.025159*(t1762*t1788 + t1945)) + 0.69051*(t1957 + t1962 + t1966 + t1977) + 0.19605*(t1957 + t1962 + t1980 + t1981 + t1985 + t1989));
  p_output1[9]=0.28330859104971495*(0.69051*(-0.15232*t1936 + t2001 + t2006 - 0.025367*t2012) + 0.19605*(t1804*t1936 + t2001 + t2006 + t1791*t2012 - 0.312342*(t1945 + t2030) - 0.025159*(t1788*t1936 + t2034)));
  p_output1[10]=0;
  p_output1[11]=0.28330859104971495*(0.69051*(-0.025367*t1774 - 0.15232*t2012 + t2042 + t2044) + 0.19605*(t1774*t1791 + t1804*t2012 - 0.025159*(t1807 + t1788*t2012) + t2042 + t2044 - 0.312342*t2065));
  p_output1[12]=0.05554264927529662*(-0.312342*(-1.*t1774*t1784 + t2030) - 0.025159*t2065 + t2012*t2075 + t1774*t2079);
  p_output1[13]=0;
  p_output1[14]=0.05554264927529662*(-0.025159*t1822 - 0.312342*(-1.*t1762*t1784 + t2063) + t1774*t2075 + t1762*t2079);
  p_output1[15]=0.28330859104971495*(0.69051*(t1864 + t1869 + t2113 + t2118) + 0.19605*(t1888 + t1895 + t1901 + t1916 + t2113 + t2118));
  p_output1[16]=0;
  p_output1[17]=0.28330859104971495*(0.69051*(t1966 + t1977 + t2129 + t2130) + 0.19605*(t1980 + t1981 + t1985 + t1989 + t2129 + t2130));
  p_output1[18]=0.05554264927529662*(-0.025157*(t1904 - 1.*t1875*t2140) - 0.312334*(t1896 - 1.*t1883*t2140) + t2140*t2144 + t1867*t2147);
  p_output1[19]=0;
  p_output1[20]=0.05554264927529662*(-0.312334*t1911 - 0.025157*(-1.*t1867*t1875 + t1983) + t1867*t2144 + t1863*t2147);
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "J_COMPosition.hh"

namespace SymFunction
{

void J_COMPosition_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
