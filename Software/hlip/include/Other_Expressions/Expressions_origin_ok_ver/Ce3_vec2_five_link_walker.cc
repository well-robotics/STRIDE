/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:37 GMT-05:00
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
  double t1716;
  double t1709;
  double t1714;
  double t1730;
  double t1756;
  double t1680;
  double t1715;
  double t1749;
  double t1754;
  double t1755;
  double t1767;
  double t1771;
  double t1782;
  double t1786;
  double t1787;
  double t1803;
  double t1804;
  double t1805;
  double t1806;
  double t1807;
  double t1874;
  double t1871;
  double t1872;
  double t1875;
  double t1873;
  double t1884;
  double t1885;
  double t1888;
  double t1889;
  double t1896;
  double t1897;
  double t1902;
  double t1905;
  double t1914;
  double t1915;
  double t1916;
  double t1917;
  double t1918;
  double t1940;
  double t1941;
  double t1942;
  double t1799;
  double t1800;
  double t1801;
  double t1830;
  double t1826;
  double t1827;
  double t1828;
  double t1829;
  double t1831;
  double t1975;
  double t1978;
  double t1979;
  double t1907;
  double t1909;
  double t1911;
  double t1931;
  double t1927;
  double t1928;
  double t1929;
  double t1930;
  double t1932;
  double t1943;
  double t1944;
  double t1945;
  double t1962;
  double t1961;
  double t1969;
  double t1948;
  double t1957;
  double t1980;
  double t1982;
  double t1984;
  double t1992;
  double t1991;
  double t1993;
  double t1988;
  double t1989;
  double t2050;
  double t2051;
  double t2052;
  double t2054;
  double t2055;
  double t2056;
  double t2070;
  double t2071;
  double t2072;
  double t2074;
  double t2075;
  double t2076;
  double t1802;
  double t1823;
  double t1824;
  double t1825;
  double t1808;
  double t1819;
  double t1820;
  double t1821;
  double t2088;
  double t2089;
  double t2090;
  double t2091;
  double t2092;
  double t1946;
  double t1947;
  double t1997;
  double t1998;
  double t1999;
  double t2000;
  double t2001;
  double t2002;
  double t2003;
  double t2004;
  double t2005;
  double t2008;
  double t2011;
  double t2017;
  double t2018;
  double t2019;
  double t2044;
  double t2045;
  double t2046;
  double t2047;
  double t2048;
  double t2049;
  double t2053;
  double t2057;
  double t2058;
  double t2060;
  double t2061;
  double t2062;
  double t2111;
  double t2112;
  double t2113;
  double t2093;
  double t2094;
  double t2095;
  double t2098;
  double t2099;
  double t2102;
  double t2103;
  double t2104;
  double t2105;
  double t2106;
  double t2107;
  double t2110;
  double t2115;
  double t2116;
  double t2120;
  double t2145;
  double t2146;
  double t2122;
  double t2148;
  double t2149;
  double t2124;
  double t1913;
  double t1924;
  double t1925;
  double t1926;
  double t1919;
  double t1920;
  double t1921;
  double t1922;
  double t2161;
  double t2162;
  double t2163;
  double t2164;
  double t2165;
  double t1986;
  double t1987;
  double t2021;
  double t2022;
  double t2023;
  double t2024;
  double t2025;
  double t2026;
  double t2027;
  double t2028;
  double t2029;
  double t2030;
  double t2031;
  double t2037;
  double t2038;
  double t2039;
  double t2064;
  double t2065;
  double t2066;
  double t2067;
  double t2068;
  double t2069;
  double t2073;
  double t2077;
  double t2078;
  double t2080;
  double t2081;
  double t2082;
  double t2184;
  double t2185;
  double t2186;
  double t2166;
  double t2167;
  double t2168;
  double t2171;
  double t2172;
  double t2175;
  double t2176;
  double t2177;
  double t2178;
  double t2179;
  double t2180;
  double t2183;
  double t2188;
  double t2189;
  double t2193;
  double t2218;
  double t2219;
  double t2195;
  double t2221;
  double t2222;
  double t2197;
  t1716 = Cos(var1[3]);
  t1709 = Cos(var1[4]);
  t1714 = Sin(var1[3]);
  t1730 = Sin(var1[4]);
  t1756 = Cos(var1[2]);
  t1680 = Sin(var1[2]);
  t1715 = -1.*t1709*t1714;
  t1749 = -1.*t1716*t1730;
  t1754 = t1715 + t1749;
  t1755 = -1.*t1680*t1754;
  t1767 = t1716*t1709;
  t1771 = -1.*t1714*t1730;
  t1782 = t1767 + t1771;
  t1786 = -1.*t1756*t1782;
  t1787 = t1755 + t1786;
  t1803 = -1.*t1709;
  t1804 = 1. + t1803;
  t1805 = 0.5*t1804;
  t1806 = 0.671885*t1709;
  t1807 = t1805 + t1806;
  t1874 = Cos(var1[5]);
  t1871 = Cos(var1[6]);
  t1872 = Sin(var1[5]);
  t1875 = Sin(var1[6]);
  t1873 = -1.*t1871*t1872;
  t1884 = -1.*t1874*t1875;
  t1885 = t1873 + t1884;
  t1888 = -1.*t1680*t1885;
  t1889 = t1874*t1871;
  t1896 = -1.*t1872*t1875;
  t1897 = t1889 + t1896;
  t1902 = -1.*t1756*t1897;
  t1905 = t1888 + t1902;
  t1914 = -1.*t1871;
  t1915 = 1. + t1914;
  t1916 = 0.5*t1915;
  t1917 = 0.671885*t1871;
  t1918 = t1916 + t1917;
  t1940 = -1.*t1716*t1680;
  t1941 = -1.*t1756*t1714;
  t1942 = t1940 + t1941;
  t1799 = -1.*t1756*t1716;
  t1800 = t1680*t1714;
  t1801 = t1799 + t1800;
  t1830 = -1.*t1680*t1782;
  t1826 = t1709*t1714;
  t1827 = t1716*t1730;
  t1828 = t1826 + t1827;
  t1829 = -1.*t1756*t1828;
  t1831 = t1829 + t1830;
  t1975 = -1.*t1874*t1680;
  t1978 = -1.*t1756*t1872;
  t1979 = t1975 + t1978;
  t1907 = -1.*t1756*t1874;
  t1909 = t1680*t1872;
  t1911 = t1907 + t1909;
  t1931 = -1.*t1680*t1897;
  t1927 = t1871*t1872;
  t1928 = t1874*t1875;
  t1929 = t1927 + t1928;
  t1930 = -1.*t1756*t1929;
  t1932 = t1930 + t1931;
  t1943 = t1756*t1716;
  t1944 = -1.*t1680*t1714;
  t1945 = t1943 + t1944;
  t1962 = t1756*t1782;
  t1961 = -1.*t1680*t1828;
  t1969 = t1961 + t1962;
  t1948 = t1756*t1754;
  t1957 = t1948 + t1830;
  t1980 = t1756*t1874;
  t1982 = -1.*t1680*t1872;
  t1984 = t1980 + t1982;
  t1992 = t1756*t1897;
  t1991 = -1.*t1680*t1929;
  t1993 = t1991 + t1992;
  t1988 = t1756*t1885;
  t1989 = t1988 + t1931;
  t2050 = t1807*t1714;
  t2051 = 0.171885*t1716*t1730;
  t2052 = t2050 + t2051;
  t2054 = t1716*t1807;
  t2055 = -0.171885*t1714*t1730;
  t2056 = t2054 + t2055;
  t2070 = t1918*t1872;
  t2071 = 0.171885*t1874*t1875;
  t2072 = t2070 + t2071;
  t2074 = t1874*t1918;
  t2075 = -0.171885*t1872*t1875;
  t2076 = t2074 + t2075;
  t1802 = -0.51185934*t1801;
  t1823 = t1807*t1730;
  t1824 = -0.171885*t1709*t1730;
  t1825 = t1823 + t1824;
  t1808 = t1807*t1709;
  t1819 = Power(t1730,2);
  t1820 = 0.171885*t1819;
  t1821 = t1808 + t1820;
  t2088 = -1.*t1716*t1709;
  t2089 = t1714*t1730;
  t2090 = t2088 + t2089;
  t2091 = t1756*t2090;
  t2092 = t1755 + t2091;
  t1946 = -6.8522*t1942*t1945;
  t1947 = -6.8522*t1942*t1801;
  t1997 = Power(t1942,2);
  t1998 = -3.4261*t1997;
  t1999 = t1716*t1680;
  t2000 = t1756*t1714;
  t2001 = t1999 + t2000;
  t2002 = -3.4261*t1942*t2001;
  t2003 = Power(t1945,2);
  t2004 = -3.4261*t2003;
  t2005 = -3.4261*t1945*t1801;
  t2008 = t1680*t1754;
  t2011 = t2008 + t1962;
  t2017 = t1756*t1828;
  t2018 = t1680*t1782;
  t2019 = t2017 + t2018;
  t2044 = Power(t1716,2);
  t2045 = 0.1494*t2044;
  t2046 = Power(t1714,2);
  t2047 = 0.1494*t2046;
  t2048 = t2045 + t2047;
  t2049 = -3.4261*t1801*t2048;
  t2053 = -1.*t2052*t1782;
  t2057 = -1.*t1754*t2056;
  t2058 = t2053 + t2057;
  t2060 = t2052*t1828;
  t2061 = t1782*t2056;
  t2062 = t2060 + t2061;
  t2111 = -1.*t1807*t1714;
  t2112 = -0.171885*t1716*t1730;
  t2113 = t2111 + t2112;
  t2093 = 0.0732367608*var2[4]*t2092;
  t2094 = -0.85216*t1825*t1957;
  t2095 = -0.85216*t1821*t2092;
  t2098 = -1.70432*t1969*t1957;
  t2099 = -1.70432*t1957*t2092;
  t2102 = -0.85216*t2011*t1969;
  t2103 = -0.85216*t1957*t2019;
  t2104 = -0.85216*t2011*t2092;
  t2105 = t1680*t2090;
  t2106 = t1948 + t2105;
  t2107 = -0.85216*t1957*t2106;
  t2110 = -0.85216*t1957*t2058;
  t2115 = t2052*t1782;
  t2116 = t1754*t2056;
  t2120 = -0.85216*t2062*t2092;
  t2145 = -0.171885*t1709*t1714;
  t2146 = t2145 + t2112;
  t2122 = -1.*t1754*t2052;
  t2148 = 0.171885*t1716*t1709;
  t2149 = t2148 + t2055;
  t2124 = -1.*t2056*t2090;
  t1913 = -0.51185934*t1911;
  t1924 = t1918*t1875;
  t1925 = -0.171885*t1871*t1875;
  t1926 = t1924 + t1925;
  t1919 = t1918*t1871;
  t1920 = Power(t1875,2);
  t1921 = 0.171885*t1920;
  t1922 = t1919 + t1921;
  t2161 = -1.*t1874*t1871;
  t2162 = t1872*t1875;
  t2163 = t2161 + t2162;
  t2164 = t1756*t2163;
  t2165 = t1888 + t2164;
  t1986 = -6.8522*t1979*t1984;
  t1987 = -6.8522*t1979*t1911;
  t2021 = Power(t1979,2);
  t2022 = -3.4261*t2021;
  t2023 = t1874*t1680;
  t2024 = t1756*t1872;
  t2025 = t2023 + t2024;
  t2026 = -3.4261*t1979*t2025;
  t2027 = Power(t1984,2);
  t2028 = -3.4261*t2027;
  t2029 = -3.4261*t1984*t1911;
  t2030 = t1680*t1885;
  t2031 = t2030 + t1992;
  t2037 = t1756*t1929;
  t2038 = t1680*t1897;
  t2039 = t2037 + t2038;
  t2064 = Power(t1874,2);
  t2065 = 0.1494*t2064;
  t2066 = Power(t1872,2);
  t2067 = 0.1494*t2066;
  t2068 = t2065 + t2067;
  t2069 = -3.4261*t1911*t2068;
  t2073 = -1.*t2072*t1897;
  t2077 = -1.*t1885*t2076;
  t2078 = t2073 + t2077;
  t2080 = t2072*t1929;
  t2081 = t1897*t2076;
  t2082 = t2080 + t2081;
  t2184 = -1.*t1918*t1872;
  t2185 = -0.171885*t1874*t1875;
  t2186 = t2184 + t2185;
  t2166 = 0.0732367608*var2[6]*t2165;
  t2167 = -0.85216*t1926*t1989;
  t2168 = -0.85216*t1922*t2165;
  t2171 = -1.70432*t1993*t1989;
  t2172 = -1.70432*t1989*t2165;
  t2175 = -0.85216*t2031*t1993;
  t2176 = -0.85216*t1989*t2039;
  t2177 = -0.85216*t2031*t2165;
  t2178 = t1680*t2163;
  t2179 = t1988 + t2178;
  t2180 = -0.85216*t1989*t2179;
  t2183 = -0.85216*t1989*t2078;
  t2188 = t2072*t1897;
  t2189 = t1885*t2076;
  t2193 = -0.85216*t2082*t2165;
  t2218 = -0.171885*t1871*t1872;
  t2219 = t2218 + t2185;
  t2195 = -1.*t1885*t2072;
  t2221 = 0.171885*t1874*t1871;
  t2222 = t2221 + t2075;
  t2197 = -1.*t2076*t2163;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(-0.5*(-0.85216*Power(t1957,2) - 0.85216*Power(t1969,2) - 0.85216*Power(t1989,2) - 0.85216*Power(t1993,2) + t1998 + t2002 + t2004 + t2005 - 0.85216*t1787*t2011 - 0.85216*t1831*t2019 + t2022 + t2026 + t2028 + t2029 - 0.85216*t1905*t2031 - 0.85216*t1932*t2039)*var2[0] - 0.5*(t1946 + t1947 - 1.70432*t1787*t1957 - 1.70432*t1831*t1969 + t1986 + t1987 - 1.70432*t1905*t1989 - 1.70432*t1932*t1993)*var2[1] - 0.5*(3.70591*t1756 + t2049 - 0.85216*t1831*t2058 - 0.85216*t1787*t2062 + t2069 - 0.85216*t1932*t2078 - 0.85216*t1905*t2082)*var2[2] - 0.5*(t1802 - 0.85216*t1787*t1821 - 0.85216*t1825*t1831)*var2[3] + 0.0732367608*t1787*var2[4] - 0.5*(t1913 - 0.85216*t1905*t1922 - 0.85216*t1926*t1932)*var2[5] + 0.0732367608*t1905*var2[6]);
  p_output1[3]=var2[1]*(t2093 - 0.5*(t1998 + t2002 + t2004 + t2005 + t2102 + t2103 + t2104 + t2107)*var2[0] - 0.5*(t1946 + t1947 + t2098 + t2099)*var2[1] - 0.5*(t2049 + t2110 - 0.85216*t1957*(t1828*t2056 + t1782*t2113 + t2115 + t2116) + t2120 - 0.85216*t1969*(-1.*t1782*t2056 - 1.*t1754*t2113 + t2122 + t2124))*var2[2] - 0.5*(t1802 + t2094 + t2095)*var2[3]);
  p_output1[4]=var2[1]*(t2093 - 0.5*(t2102 + t2103 + t2104 + t2107)*var2[0] - 0.5*(t2098 + t2099)*var2[1] - 0.5*(t2110 + t2120 - 0.85216*t1969*(t2122 + t2124 - 1.*t1754*t2146 - 1.*t1782*t2149) - 0.85216*t1957*(t2115 + t2116 + t1782*t2146 + t1828*t2149))*var2[2] - 0.5*(-0.85216*(0.171885*t1709*t1730 - 1.*t1730*t1807)*t1957 - 0.85216*(-0.171885*Power(t1709,2) + t1808)*t1969 + t2094 + t2095)*var2[3]);
  p_output1[5]=var2[1]*(t2166 - 0.5*(t2022 + t2026 + t2028 + t2029 + t2175 + t2176 + t2177 + t2180)*var2[0] - 0.5*(t1986 + t1987 + t2171 + t2172)*var2[1] - 0.5*(t2069 + t2183 - 0.85216*t1989*(t1929*t2076 + t1897*t2186 + t2188 + t2189) + t2193 - 0.85216*t1993*(-1.*t1897*t2076 - 1.*t1885*t2186 + t2195 + t2197))*var2[2] - 0.5*(t1913 + t2167 + t2168)*var2[5]);
  p_output1[6]=var2[1]*(t2166 - 0.5*(t2175 + t2176 + t2177 + t2180)*var2[0] - 0.5*(t2171 + t2172)*var2[1] - 0.5*(t2183 + t2193 - 0.85216*t1993*(t2195 + t2197 - 1.*t1885*t2219 - 1.*t1897*t2222) - 0.85216*t1989*(t2188 + t2189 + t1897*t2219 + t1929*t2222))*var2[2] - 0.5*(-0.85216*(0.171885*t1871*t1875 - 1.*t1875*t1918)*t1989 - 0.85216*(-0.171885*Power(t1871,2) + t1919)*t1993 + t2167 + t2168)*var2[5]);
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

#include "Ce3_vec2_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
