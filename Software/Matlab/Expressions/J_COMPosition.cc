/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:22 GMT-05:00
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
  double t1699;
  double t1641;
  double t1765;
  double t1768;
  double t1773;
  double t1782;
  double t1779;
  double t1783;
  double t1791;
  double t1795;
  double t1796;
  double t1798;
  double t1801;
  double t1802;
  double t1814;
  double t1818;
  double t1819;
  double t1826;
  double t1827;
  double t1831;
  double t1841;
  double t1842;
  double t1843;
  double t1845;
  double t1902;
  double t1903;
  double t1904;
  double t1907;
  double t1905;
  double t1908;
  double t1909;
  double t1912;
  double t1913;
  double t1914;
  double t1915;
  double t1917;
  double t1922;
  double t1923;
  double t1924;
  double t1918;
  double t1919;
  double t1920;
  double t1932;
  double t1934;
  double t1937;
  double t1940;
  double t1987;
  double t1988;
  double t1990;
  double t1991;
  double t1993;
  double t1844;
  double t1846;
  double t1848;
  double t1858;
  double t1860;
  double t1861;
  double t1884;
  double t2014;
  double t2019;
  double t1939;
  double t1941;
  double t1942;
  double t2026;
  double t2027;
  double t2033;
  double t1947;
  double t1950;
  double t1951;
  double t1959;
  double t2056;
  double t2057;
  double t2058;
  double t2059;
  double t2060;
  double t2061;
  double t2063;
  double t2064;
  double t2065;
  double t2069;
  double t2002;
  double t2099;
  double t2101;
  double t1864;
  double t2091;
  double t2087;
  double t2120;
  double t2122;
  double t2131;
  double t2132;
  double t2134;
  double t2135;
  double t2136;
  double t1870;
  double t1879;
  double t1921;
  double t1926;
  double t2164;
  double t2165;
  double t2166;
  double t2170;
  double t2172;
  double t2173;
  double t2175;
  double t1945;
  double t1952;
  double t1953;
  double t1954;
  double t1955;
  double t1958;
  double t1961;
  double t1968;
  double t1973;
  double t2023;
  double t2034;
  double t2186;
  double t2187;
  double t2037;
  double t2038;
  double t2039;
  double t2040;
  double t2041;
  double t2042;
  double t2044;
  double t2045;
  double t2046;
  double t2195;
  double t2196;
  double t2197;
  double t2199;
  double t2200;
  double t2201;
  double t2203;
  double t2204;
  t1699 = Sin(var1[2]);
  t1641 = Cos(var1[2]);
  t1765 = Cos(var1[3]);
  t1768 = -1.*t1765;
  t1773 = 1. + t1768;
  t1782 = Sin(var1[3]);
  t1779 = -0.0265*t1773;
  t1783 = -0.0695*t1782;
  t1791 = t1779 + t1783;
  t1795 = -1.*t1699*t1791;
  t1796 = -0.0695*t1773;
  t1798 = 0.0265*t1782;
  t1801 = t1796 + t1798;
  t1802 = t1641*t1801;
  t1814 = -1.*t1765*t1699;
  t1818 = t1641*t1782;
  t1819 = t1814 + t1818;
  t1826 = t1641*t1765;
  t1827 = t1699*t1782;
  t1831 = t1826 + t1827;
  t1841 = Cos(var1[4]);
  t1842 = -1.*t1841;
  t1843 = 1. + t1842;
  t1845 = Sin(var1[4]);
  t1902 = Cos(var1[5]);
  t1903 = -1.*t1902;
  t1904 = 1. + t1903;
  t1907 = Sin(var1[5]);
  t1905 = -0.0695*t1904;
  t1908 = -0.0265*t1907;
  t1909 = t1905 + t1908;
  t1912 = t1641*t1909;
  t1913 = -0.0265*t1904;
  t1914 = 0.0695*t1907;
  t1915 = t1913 + t1914;
  t1917 = -1.*t1699*t1915;
  t1922 = t1641*t1902;
  t1923 = -1.*t1699*t1907;
  t1924 = t1922 + t1923;
  t1918 = -1.*t1902*t1699;
  t1919 = -1.*t1641*t1907;
  t1920 = t1918 + t1919;
  t1932 = Cos(var1[6]);
  t1934 = -1.*t1932;
  t1937 = 1. + t1934;
  t1940 = Sin(var1[6]);
  t1987 = -1.*t1641*t1791;
  t1988 = -1.*t1699*t1801;
  t1990 = -1.*t1641*t1765;
  t1991 = -1.*t1699*t1782;
  t1993 = t1990 + t1991;
  t1844 = -0.0265*t1843;
  t1846 = -0.2375*t1845;
  t1848 = t1844 + t1846;
  t1858 = -0.2375*t1843;
  t1860 = 0.0265*t1845;
  t1861 = t1858 + t1860;
  t1884 = t1841*t1819;
  t2014 = -1.*t1699*t1909;
  t2019 = -1.*t1641*t1915;
  t1939 = -0.2375*t1937;
  t1941 = -0.0265*t1940;
  t1942 = t1939 + t1941;
  t2026 = -1.*t1641*t1902;
  t2027 = t1699*t1907;
  t2033 = t2026 + t2027;
  t1947 = -0.0265*t1937;
  t1950 = 0.2375*t1940;
  t1951 = t1947 + t1950;
  t1959 = t1932*t1920;
  t2056 = 0.0265*t1765;
  t2057 = t2056 + t1783;
  t2058 = t1699*t2057;
  t2059 = -0.0695*t1765;
  t2060 = -0.0265*t1782;
  t2061 = t2059 + t2060;
  t2063 = t1641*t2061;
  t2064 = t1765*t1699;
  t2065 = -1.*t1641*t1782;
  t2069 = t2064 + t2065;
  t2002 = t1841*t1993;
  t2099 = t1641*t2057;
  t2101 = -1.*t1699*t2061;
  t1864 = t1841*t1831;
  t2091 = t1841*t2069;
  t2087 = -1.*t2069*t1845;
  t2120 = -1.*t1831*t1845;
  t2122 = t2091 + t2120;
  t2131 = 0.0265*t1841;
  t2132 = t2131 + t1846;
  t2134 = -0.2375*t1841;
  t2135 = -0.0265*t1845;
  t2136 = t2134 + t2135;
  t1870 = -1.*t1819*t1845;
  t1879 = t1864 + t1870;
  t1921 = -0.025413*t1920;
  t1926 = -0.15232*t1924;
  t2164 = -0.0265*t1902;
  t2165 = -0.0695*t1907;
  t2166 = t2164 + t2165;
  t2170 = t1699*t2166;
  t2172 = 0.0695*t1902;
  t2173 = t2172 + t1908;
  t2175 = t1641*t2173;
  t1945 = t1924*t1942;
  t1952 = t1920*t1951;
  t1953 = t1932*t1924;
  t1954 = t1920*t1940;
  t1955 = t1953 + t1954;
  t1958 = -0.314506*t1955;
  t1961 = -1.*t1924*t1940;
  t1968 = t1959 + t1961;
  t1973 = -0.025226*t1968;
  t2023 = -0.15232*t1920;
  t2034 = -0.025413*t2033;
  t2186 = t1641*t2166;
  t2187 = -1.*t1699*t2173;
  t2037 = t1920*t1942;
  t2038 = t2033*t1951;
  t2039 = t1932*t2033;
  t2040 = -1.*t1920*t1940;
  t2041 = t2039 + t2040;
  t2042 = -0.025226*t2041;
  t2044 = t2033*t1940;
  t2045 = t1959 + t2044;
  t2046 = -0.314506*t2045;
  t2195 = t1902*t1699;
  t2196 = t1641*t1907;
  t2197 = t2195 + t2196;
  t2199 = -0.0265*t1932;
  t2200 = -0.2375*t1940;
  t2201 = t2199 + t2200;
  t2203 = 0.2375*t1932;
  t2204 = t2203 + t1941;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=1.;
  p_output1[6]=0.2996793431028799*(1.5566*(-0.046589*t1641 + 0.026461*t1699) + 0.69051*(t1795 + t1802 - 0.025367*t1819 - 0.15232*t1831) + 0.19964*(t1795 + t1802 + t1819*t1848 + t1831*t1861 - 0.314514*t1879 - 0.025229*(t1831*t1845 + t1884)) + 0.69051*(t1912 + t1917 + t1921 + t1926) + 0.19964*(t1912 + t1917 + t1945 + t1952 + t1958 + t1973));
  p_output1[7]=0;
  p_output1[8]=0.2996793431028799*(1.5566*(0.026461*t1641 + 0.046589*t1699) + 0.69051*(-0.15232*t1819 + t1987 + t1988 - 0.025367*t1993) + 0.19964*(t1819*t1861 + t1987 + t1988 + t1848*t1993 - 0.314514*(t1884 - 1.*t1845*t1993) - 0.025229*(t1819*t1845 + t2002)) + 0.69051*(t2014 + t2019 + t2023 + t2034) + 0.19964*(t2014 + t2019 + t2037 + t2038 + t2042 + t2046));
  p_output1[9]=0.2996793431028799*(0.69051*(-0.15232*t1993 + t2058 + t2063 - 0.025367*t2069) + 0.19964*(t1861*t1993 + t2058 + t2063 + t1848*t2069 - 0.314514*(t2002 + t2087) - 0.025229*(t1845*t1993 + t2091)));
  p_output1[10]=0;
  p_output1[11]=0.2996793431028799*(0.69051*(-0.025367*t1831 - 0.15232*t2069 + t2099 + t2101) + 0.19964*(t1831*t1848 + t1861*t2069 - 0.025229*(t1864 + t1845*t2069) + t2099 + t2101 - 0.314514*t2122));
  p_output1[12]=0.05982798405705895*(-0.314514*(-1.*t1831*t1841 + t2087) - 0.025229*t2122 + t2069*t2132 + t1831*t2136);
  p_output1[13]=0;
  p_output1[14]=0.05982798405705895*(-0.025229*t1879 - 0.314514*(-1.*t1819*t1841 + t2120) + t1831*t2132 + t1819*t2136);
  p_output1[15]=0.2996793431028799*(0.69051*(t1921 + t1926 + t2170 + t2175) + 0.19964*(t1945 + t1952 + t1958 + t1973 + t2170 + t2175));
  p_output1[16]=0;
  p_output1[17]=0.2996793431028799*(0.69051*(t2023 + t2034 + t2186 + t2187) + 0.19964*(t2037 + t2038 + t2042 + t2046 + t2186 + t2187));
  p_output1[18]=0.05982798405705895*(-0.025226*(t1961 - 1.*t1932*t2197) - 0.314506*(t1953 - 1.*t1940*t2197) + t2197*t2201 + t1924*t2204);
  p_output1[19]=0;
  p_output1[20]=0.05982798405705895*(-0.314506*t1968 - 0.025226*(-1.*t1924*t1932 + t2040) + t1924*t2201 + t1920*t2204);
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
