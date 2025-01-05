/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:25 GMT-05:00
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
  double t1750;
  double t1679;
  double t1700;
  double t1779;
  double t1798;
  double t1818;
  double t1832;
  double t1836;
  double t1837;
  double t1746;
  double t1783;
  double t1796;
  double t1893;
  double t1899;
  double t1905;
  double t1908;
  double t1912;
  double t1913;
  double t1914;
  double t1921;
  double t1825;
  double t1926;
  double t1927;
  double t1928;
  double t1945;
  double t1947;
  double t1860;
  double t1864;
  double t1989;
  double t1995;
  double t2001;
  double t2004;
  double t2009;
  double t2010;
  double t2011;
  double t1994;
  double t1999;
  double t2000;
  double t2046;
  double t2047;
  double t2048;
  double t2049;
  double t2052;
  double t2056;
  double t2057;
  double t2058;
  double t2023;
  double t2059;
  double t2060;
  double t2063;
  double t2084;
  double t2085;
  double t2034;
  double t2039;
  double t2122;
  double t2123;
  double t2124;
  double t2126;
  double t2127;
  double t2128;
  double t2131;
  double t2132;
  double t2133;
  double t2134;
  double t2135;
  double t2136;
  double t1950;
  double t1951;
  double t1954;
  double t1955;
  double t1958;
  double t2165;
  double t2166;
  double t2064;
  double t2075;
  double t2170;
  double t2172;
  double t2175;
  double t2176;
  double t2179;
  double t2180;
  double t2181;
  double t2086;
  double t2088;
  double t2089;
  double t2091;
  double t2092;
  double t2093;
  double t2094;
  double t2095;
  double t2096;
  double t2097;
  double t2098;
  double t2099;
  double t2102;
  double t1802;
  double t1826;
  double t1844;
  double t1846;
  double t1848;
  double t2146;
  double t2147;
  double t2149;
  double t1879;
  double t2203;
  double t2204;
  double t1882;
  double t2002;
  double t2005;
  double t2006;
  double t2013;
  double t2027;
  double t2240;
  double t2241;
  double t2245;
  double t2247;
  double t2248;
  double t2229;
  double t2281;
  double t2282;
  double t1968;
  double t2151;
  double t2243;
  double t2249;
  double t2305;
  double t2307;
  double t2257;
  double t2260;
  double t2261;
  double t2262;
  double t2263;
  double t2264;
  double t2265;
  double t2266;
  t1750 = Cos(var1[2]);
  t1679 = Cos(var1[3]);
  t1700 = Sin(var1[2]);
  t1779 = Sin(var1[3]);
  t1798 = Cos(var1[4]);
  t1818 = Sin(var1[4]);
  t1832 = t1750*t1679;
  t1836 = t1700*t1779;
  t1837 = t1832 + t1836;
  t1746 = t1679*t1700;
  t1783 = -1.*t1750*t1779;
  t1796 = t1746 + t1783;
  t1893 = 0.0265*t1679;
  t1899 = -0.0695*t1779;
  t1905 = t1893 + t1899;
  t1908 = t1700*t1905;
  t1912 = -0.0695*t1679;
  t1913 = -0.0265*t1779;
  t1914 = t1912 + t1913;
  t1921 = t1750*t1914;
  t1825 = -0.2375*t1818;
  t1926 = -1.*t1750*t1679;
  t1927 = -1.*t1700*t1779;
  t1928 = t1926 + t1927;
  t1945 = -1.*t1798;
  t1947 = 1. + t1945;
  t1860 = -1.*t1796*t1818;
  t1864 = t1798*t1796;
  t1989 = Cos(var1[5]);
  t1995 = Sin(var1[5]);
  t2001 = Cos(var1[6]);
  t2004 = Sin(var1[6]);
  t2009 = t1750*t1989;
  t2010 = -1.*t1700*t1995;
  t2011 = t2009 + t2010;
  t1994 = t1989*t1700;
  t1999 = t1750*t1995;
  t2000 = t1994 + t1999;
  t2046 = -0.0265*t1989;
  t2047 = -0.0695*t1995;
  t2048 = t2046 + t2047;
  t2049 = t1700*t2048;
  t2052 = 0.0695*t1989;
  t2056 = -0.0265*t1995;
  t2057 = t2052 + t2056;
  t2058 = t1750*t2057;
  t2023 = -0.0265*t2004;
  t2059 = -1.*t1989*t1700;
  t2060 = -1.*t1750*t1995;
  t2063 = t2059 + t2060;
  t2084 = -1.*t2001;
  t2085 = 1. + t2084;
  t2034 = t2001*t2011;
  t2039 = -1.*t2011*t2004;
  t2122 = -1.*t1679;
  t2123 = 1. + t2122;
  t2124 = -0.0265*t2123;
  t2126 = t2124 + t1899;
  t2127 = -1.*t1700*t2126;
  t2128 = -0.0695*t2123;
  t2131 = 0.0265*t1779;
  t2132 = t2128 + t2131;
  t2133 = t1750*t2132;
  t2134 = -1.*t1679*t1700;
  t2135 = t1750*t1779;
  t2136 = t2134 + t2135;
  t1950 = -0.0265*t1947;
  t1951 = t1950 + t1825;
  t1954 = -0.2375*t1947;
  t1955 = 0.0265*t1818;
  t1958 = t1954 + t1955;
  t2165 = -1.*t1989;
  t2166 = 1. + t2165;
  t2064 = -0.025413*t2063;
  t2075 = -0.15232*t2011;
  t2170 = -0.0695*t2166;
  t2172 = t2170 + t2056;
  t2175 = t1750*t2172;
  t2176 = -0.0265*t2166;
  t2179 = 0.0695*t1995;
  t2180 = t2176 + t2179;
  t2181 = -1.*t1700*t2180;
  t2086 = -0.2375*t2085;
  t2088 = t2086 + t2023;
  t2089 = t2011*t2088;
  t2091 = -0.0265*t2085;
  t2092 = 0.2375*t2004;
  t2093 = t2091 + t2092;
  t2094 = t2063*t2093;
  t2095 = t2063*t2004;
  t2096 = t2034 + t2095;
  t2097 = -0.314506*t2096;
  t2098 = t2001*t2063;
  t2099 = t2098 + t2039;
  t2102 = -0.025226*t2099;
  t1802 = 0.0265*t1798;
  t1826 = t1802 + t1825;
  t1844 = -0.2375*t1798;
  t1846 = -0.0265*t1818;
  t1848 = t1844 + t1846;
  t2146 = t1798*t1837;
  t2147 = -1.*t2136*t1818;
  t2149 = t2146 + t2147;
  t1879 = -1.*t1837*t1818;
  t2203 = t1750*t1905;
  t2204 = -1.*t1700*t1914;
  t1882 = t1864 + t1879;
  t2002 = -0.0265*t2001;
  t2005 = -0.2375*t2004;
  t2006 = t2002 + t2005;
  t2013 = 0.2375*t2001;
  t2027 = t2013 + t2023;
  t2240 = t1750*t2048;
  t2241 = -1.*t1700*t2057;
  t2245 = -1.*t1750*t1989;
  t2247 = t1700*t1995;
  t2248 = t2245 + t2247;
  t2229 = -1.*t2063*t2004;
  t2281 = -1.*t1750*t2126;
  t2282 = -1.*t1700*t2132;
  t1968 = t1798*t1928;
  t2151 = t1798*t2136;
  t2243 = -0.15232*t2063;
  t2249 = -0.025413*t2248;
  t2305 = -1.*t1700*t2172;
  t2307 = -1.*t1750*t2180;
  t2257 = t2063*t2088;
  t2260 = t2248*t2093;
  t2261 = t2001*t2248;
  t2262 = t2261 + t2229;
  t2263 = -0.025226*t2262;
  t2264 = t2248*t2004;
  t2265 = t2098 + t2264;
  t2266 = -0.314506*t2265;
  p_output1[0]=var2[0] + 0.2996793431028799*(1.5566*(0.026461*t1700 - 0.046589*t1750) + 0.69051*(-0.15232*t1837 + t2127 + t2133 - 0.025367*t2136) + 0.19964*(t1837*t1958 + t2127 + t2133 + t1951*t2136 - 0.314514*t2149 - 0.025229*(t1818*t1837 + t2151)) + 0.69051*(t2064 + t2075 + t2175 + t2181) + 0.19964*(t2089 + t2094 + t2097 + t2102 + t2175 + t2181))*var2[2] + 0.2996793431028799*(0.69051*(-0.025367*t1796 + t1908 + t1921 - 0.15232*t1928) + 0.19964*(t1908 + t1921 - 0.025229*(t1864 + t1818*t1928) + t1796*t1951 + t1928*t1958 - 0.314514*(t1860 + t1968)))*var2[3] + 0.05982798405705895*(t1796*t1826 + t1837*t1848 - 0.314514*(-1.*t1798*t1837 + t1860) - 0.025229*t1882)*var2[4] + 0.2996793431028799*(0.69051*(t2049 + t2058 + t2064 + t2075) + 0.19964*(t2049 + t2058 + t2089 + t2094 + t2097 + t2102))*var2[5] + 0.05982798405705895*(t2000*t2006 + t2011*t2027 - 0.314506*(-1.*t2000*t2004 + t2034) - 0.025226*(-1.*t2000*t2001 + t2039))*var2[6];
  p_output1[1]=0;
  p_output1[2]=var2[1] + 0.2996793431028799*(1.5566*(0.046589*t1700 + 0.026461*t1750) + 0.69051*(-0.025367*t1928 - 0.15232*t2136 + t2281 + t2282) + 0.19964*(t1928*t1951 + t1958*t2136 - 0.025229*(t1968 + t1818*t2136) - 0.314514*(-1.*t1818*t1928 + t2151) + t2281 + t2282) + 0.69051*(t2243 + t2249 + t2305 + t2307) + 0.19964*(t2257 + t2260 + t2263 + t2266 + t2305 + t2307))*var2[2] + 0.2996793431028799*(0.69051*(-0.15232*t1796 - 0.025367*t1837 + t2203 + t2204) + 0.19964*(-0.314514*t1882 + t1837*t1951 + t1796*t1958 - 0.025229*(t1796*t1818 + t2146) + t2203 + t2204))*var2[3] + 0.05982798405705895*(t1826*t1837 + t1848*t2136 - 0.314514*(t1879 - 1.*t1798*t2136) - 0.025229*t2149)*var2[4] + 0.2996793431028799*(0.69051*(t2240 + t2241 + t2243 + t2249) + 0.19964*(t2240 + t2241 + t2257 + t2260 + t2263 + t2266))*var2[5] + 0.05982798405705895*(t2006*t2011 + t2027*t2063 - 0.314506*t2099 - 0.025226*(-1.*t2001*t2011 + t2229))*var2[6];
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
