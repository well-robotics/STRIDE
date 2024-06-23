/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:38 GMT-05:00
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
  double t1722;
  double t1643;
  double t1689;
  double t1726;
  double t1745;
  double t1768;
  double t1779;
  double t1780;
  double t1787;
  double t1693;
  double t1739;
  double t1741;
  double t1842;
  double t1848;
  double t1851;
  double t1855;
  double t1856;
  double t1857;
  double t1864;
  double t1865;
  double t1769;
  double t1870;
  double t1871;
  double t1882;
  double t1890;
  double t1893;
  double t1804;
  double t1822;
  double t1937;
  double t1942;
  double t1945;
  double t1948;
  double t1953;
  double t1954;
  double t1956;
  double t1938;
  double t1943;
  double t1944;
  double t1990;
  double t1991;
  double t1992;
  double t1995;
  double t1999;
  double t2000;
  double t2001;
  double t2002;
  double t1970;
  double t2003;
  double t2006;
  double t2007;
  double t2028;
  double t2029;
  double t1978;
  double t1984;
  double t2066;
  double t2067;
  double t2069;
  double t2070;
  double t2071;
  double t2074;
  double t2075;
  double t2076;
  double t2077;
  double t2078;
  double t2079;
  double t2080;
  double t1894;
  double t1895;
  double t1898;
  double t1901;
  double t1902;
  double t2109;
  double t2113;
  double t2018;
  double t2019;
  double t2115;
  double t2118;
  double t2119;
  double t2122;
  double t2123;
  double t2124;
  double t2127;
  double t2031;
  double t2032;
  double t2034;
  double t2035;
  double t2036;
  double t2037;
  double t2038;
  double t2039;
  double t2040;
  double t2041;
  double t2042;
  double t2045;
  double t2046;
  double t1761;
  double t1770;
  double t1789;
  double t1791;
  double t1797;
  double t2090;
  double t2092;
  double t2093;
  double t1825;
  double t2148;
  double t2149;
  double t1828;
  double t1947;
  double t1949;
  double t1950;
  double t1966;
  double t1976;
  double t2186;
  double t2188;
  double t2191;
  double t2192;
  double t2197;
  double t2175;
  double t2226;
  double t2228;
  double t1916;
  double t2097;
  double t2190;
  double t2199;
  double t2251;
  double t2252;
  double t2204;
  double t2205;
  double t2206;
  double t2207;
  double t2208;
  double t2209;
  double t2211;
  double t2213;
  t1722 = Cos(var1[2]);
  t1643 = Cos(var1[3]);
  t1689 = Sin(var1[2]);
  t1726 = Sin(var1[3]);
  t1745 = Cos(var1[4]);
  t1768 = Sin(var1[4]);
  t1779 = t1722*t1643;
  t1780 = t1689*t1726;
  t1787 = t1779 + t1780;
  t1693 = t1643*t1689;
  t1739 = -1.*t1722*t1726;
  t1741 = t1693 + t1739;
  t1842 = 0.0265*t1643;
  t1848 = -0.0695*t1726;
  t1851 = t1842 + t1848;
  t1855 = t1689*t1851;
  t1856 = -0.0695*t1643;
  t1857 = -0.0265*t1726;
  t1864 = t1856 + t1857;
  t1865 = t1722*t1864;
  t1769 = -0.2375*t1768;
  t1870 = -1.*t1722*t1643;
  t1871 = -1.*t1689*t1726;
  t1882 = t1870 + t1871;
  t1890 = -1.*t1745;
  t1893 = 1. + t1890;
  t1804 = -1.*t1741*t1768;
  t1822 = t1745*t1741;
  t1937 = Cos(var1[5]);
  t1942 = Sin(var1[5]);
  t1945 = Cos(var1[6]);
  t1948 = Sin(var1[6]);
  t1953 = t1722*t1937;
  t1954 = -1.*t1689*t1942;
  t1956 = t1953 + t1954;
  t1938 = t1937*t1689;
  t1943 = t1722*t1942;
  t1944 = t1938 + t1943;
  t1990 = -0.0265*t1937;
  t1991 = -0.0695*t1942;
  t1992 = t1990 + t1991;
  t1995 = t1689*t1992;
  t1999 = 0.0695*t1937;
  t2000 = -0.0265*t1942;
  t2001 = t1999 + t2000;
  t2002 = t1722*t2001;
  t1970 = -0.0265*t1948;
  t2003 = -1.*t1937*t1689;
  t2006 = -1.*t1722*t1942;
  t2007 = t2003 + t2006;
  t2028 = -1.*t1945;
  t2029 = 1. + t2028;
  t1978 = t1945*t1956;
  t1984 = -1.*t1956*t1948;
  t2066 = -1.*t1643;
  t2067 = 1. + t2066;
  t2069 = -0.0265*t2067;
  t2070 = t2069 + t1848;
  t2071 = -1.*t1689*t2070;
  t2074 = -0.0695*t2067;
  t2075 = 0.0265*t1726;
  t2076 = t2074 + t2075;
  t2077 = t1722*t2076;
  t2078 = -1.*t1643*t1689;
  t2079 = t1722*t1726;
  t2080 = t2078 + t2079;
  t1894 = -0.0265*t1893;
  t1895 = t1894 + t1769;
  t1898 = -0.2375*t1893;
  t1901 = 0.0265*t1768;
  t1902 = t1898 + t1901;
  t2109 = -1.*t1937;
  t2113 = 1. + t2109;
  t2018 = -0.025413*t2007;
  t2019 = -0.15232*t1956;
  t2115 = -0.0695*t2113;
  t2118 = t2115 + t2000;
  t2119 = t1722*t2118;
  t2122 = -0.0265*t2113;
  t2123 = 0.0695*t1942;
  t2124 = t2122 + t2123;
  t2127 = -1.*t1689*t2124;
  t2031 = -0.2375*t2029;
  t2032 = t2031 + t1970;
  t2034 = t1956*t2032;
  t2035 = -0.0265*t2029;
  t2036 = 0.2375*t1948;
  t2037 = t2035 + t2036;
  t2038 = t2007*t2037;
  t2039 = t2007*t1948;
  t2040 = t1978 + t2039;
  t2041 = -0.312334*t2040;
  t2042 = t1945*t2007;
  t2045 = t2042 + t1984;
  t2046 = -0.025157*t2045;
  t1761 = 0.0265*t1745;
  t1770 = t1761 + t1769;
  t1789 = -0.2375*t1745;
  t1791 = -0.0265*t1768;
  t1797 = t1789 + t1791;
  t2090 = t1745*t1787;
  t2092 = -1.*t2080*t1768;
  t2093 = t2090 + t2092;
  t1825 = -1.*t1787*t1768;
  t2148 = t1722*t1851;
  t2149 = -1.*t1689*t1864;
  t1828 = t1822 + t1825;
  t1947 = -0.0265*t1945;
  t1949 = -0.2375*t1948;
  t1950 = t1947 + t1949;
  t1966 = 0.2375*t1945;
  t1976 = t1966 + t1970;
  t2186 = t1722*t1992;
  t2188 = -1.*t1689*t2001;
  t2191 = -1.*t1722*t1937;
  t2192 = t1689*t1942;
  t2197 = t2191 + t2192;
  t2175 = -1.*t2007*t1948;
  t2226 = -1.*t1722*t2070;
  t2228 = -1.*t1689*t2076;
  t1916 = t1745*t1882;
  t2097 = t1745*t2080;
  t2190 = -0.15232*t2007;
  t2199 = -0.025413*t2197;
  t2251 = -1.*t1689*t2118;
  t2252 = -1.*t1722*t2124;
  t2204 = t2007*t2032;
  t2205 = t2197*t2037;
  t2206 = t1945*t2197;
  t2207 = t2206 + t2175;
  t2208 = -0.025157*t2207;
  t2209 = t2197*t1948;
  t2211 = t2042 + t2209;
  t2213 = -0.312334*t2211;
  p_output1[0]=0.9999999999999998*var2[0] + 0.28330859104971495*(1.7566*(0.026461*t1689 - 0.046589*t1722) + 0.69051*(-0.15232*t1787 + t2071 + t2077 - 0.025367*t2080) + 0.19605*(t1787*t1902 + t2071 + t2077 + t1895*t2080 - 0.312342*t2093 - 0.025159*(t1768*t1787 + t2097)) + 0.69051*(t2018 + t2019 + t2119 + t2127) + 0.19605*(t2034 + t2038 + t2041 + t2046 + t2119 + t2127))*var2[2] + 0.28330859104971495*(0.69051*(-0.025367*t1741 + t1855 + t1865 - 0.15232*t1882) + 0.19605*(t1855 + t1865 - 0.025159*(t1822 + t1768*t1882) + t1741*t1895 + t1882*t1902 - 0.312342*(t1804 + t1916)))*var2[3] + 0.05554264927529662*(t1741*t1770 + t1787*t1797 - 0.312342*(-1.*t1745*t1787 + t1804) - 0.025159*t1828)*var2[4] + 0.28330859104971495*(0.69051*(t1995 + t2002 + t2018 + t2019) + 0.19605*(t1995 + t2002 + t2034 + t2038 + t2041 + t2046))*var2[5] + 0.05554264927529662*(t1944*t1950 + t1956*t1976 - 0.312334*(-1.*t1944*t1948 + t1978) - 0.025157*(-1.*t1944*t1945 + t1984))*var2[6];
  p_output1[1]=0;
  p_output1[2]=0.9999999999999998*var2[1] + 0.28330859104971495*(1.7566*(0.046589*t1689 + 0.026461*t1722) + 0.69051*(-0.025367*t1882 - 0.15232*t2080 + t2226 + t2228) + 0.19605*(t1882*t1895 + t1902*t2080 - 0.025159*(t1916 + t1768*t2080) - 0.312342*(-1.*t1768*t1882 + t2097) + t2226 + t2228) + 0.69051*(t2190 + t2199 + t2251 + t2252) + 0.19605*(t2204 + t2205 + t2208 + t2213 + t2251 + t2252))*var2[2] + 0.28330859104971495*(0.69051*(-0.15232*t1741 - 0.025367*t1787 + t2148 + t2149) + 0.19605*(-0.312342*t1828 + t1787*t1895 + t1741*t1902 - 0.025159*(t1741*t1768 + t2090) + t2148 + t2149))*var2[3] + 0.05554264927529662*(t1770*t1787 + t1797*t2080 - 0.312342*(t1825 - 1.*t1745*t2080) - 0.025159*t2093)*var2[4] + 0.28330859104971495*(0.69051*(t2186 + t2188 + t2190 + t2199) + 0.19605*(t2186 + t2188 + t2204 + t2205 + t2208 + t2213))*var2[5] + 0.05554264927529662*(t1950*t1956 + t1976*t2007 - 0.312334*t2045 - 0.025157*(-1.*t1945*t1956 + t2175))*var2[6];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "COM_velocity.hh"

namespace SymFunction
{

void COM_velocity_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
