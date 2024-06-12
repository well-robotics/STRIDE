/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:39 GMT-05:00
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
  double t1693;
  double t1769;
  double t1622;
  double t1770;
  double t1780;
  double t1791;
  double t1761;
  double t1775;
  double t1779;
  double t1806;
  double t1807;
  double t1822;
  double t1888;
  double t1894;
  double t1895;
  double t1897;
  double t1898;
  double t1901;
  double t1902;
  double t1911;
  double t1797;
  double t1916;
  double t1918;
  double t1919;
  double t1928;
  double t1929;
  double t1842;
  double t1857;
  double t1987;
  double t1989;
  double t1992;
  double t1999;
  double t1988;
  double t1990;
  double t1991;
  double t2019;
  double t2022;
  double t2027;
  double t2056;
  double t2057;
  double t2058;
  double t2060;
  double t2061;
  double t2062;
  double t2065;
  double t2069;
  double t2032;
  double t2075;
  double t2077;
  double t2079;
  double t2086;
  double t2088;
  double t2037;
  double t2040;
  double t2138;
  double t2139;
  double t2142;
  double t2143;
  double t2144;
  double t2145;
  double t2146;
  double t2147;
  double t2148;
  double t2154;
  double t2155;
  double t2156;
  double t1932;
  double t1938;
  double t1947;
  double t1949;
  double t1950;
  double t2180;
  double t2183;
  double t2074;
  double t2082;
  double t2184;
  double t2186;
  double t2190;
  double t2192;
  double t2197;
  double t2199;
  double t2200;
  double t2089;
  double t2090;
  double t2093;
  double t2094;
  double t2099;
  double t2101;
  double t2104;
  double t2107;
  double t2108;
  double t2115;
  double t2119;
  double t2122;
  double t2123;
  double t1789;
  double t1801;
  double t1828;
  double t1833;
  double t1834;
  double t1848;
  double t2169;
  double t2170;
  double t2172;
  double t2224;
  double t2225;
  double t1851;
  double t1995;
  double t2000;
  double t2002;
  double t2031;
  double t2034;
  double t2270;
  double t2271;
  double t2272;
  double t2273;
  double t2275;
  double t2260;
  double t2309;
  double t2311;
  double t2164;
  double t1976;
  double t2277;
  double t2279;
  double t2328;
  double t2329;
  double t2284;
  double t2285;
  double t2287;
  double t2288;
  double t2289;
  double t2293;
  double t2294;
  double t2295;
  double t2318;
  double t2319;
  double t2151;
  double t2157;
  double t2350;
  double t2351;
  double t2352;
  double t2355;
  double t2356;
  double t2161;
  double t2162;
  double t2165;
  double t2166;
  double t2167;
  double t2174;
  double t1922;
  double t1925;
  double t1926;
  double t1927;
  double t1943;
  double t1952;
  double t1953;
  double t1966;
  double t1970;
  double t1977;
  double t1979;
  double t1980;
  double t1981;
  double t1982;
  double t2312;
  double t2313;
  double t2383;
  double t2384;
  double t2316;
  double t2317;
  double t2320;
  double t2323;
  double t2324;
  double t2325;
  double t2229;
  double t2230;
  double t2231;
  double t2233;
  double t2234;
  double t2235;
  double t2236;
  double t2237;
  double t2238;
  double t2243;
  double t2246;
  double t2247;
  double t2248;
  double t2339;
  double t2341;
  double t2342;
  double t2343;
  double t2345;
  double t2346;
  double t2347;
  double t1803;
  double t1836;
  double t1855;
  double t1856;
  double t1865;
  double t1869;
  double t1884;
  double t2376;
  double t2377;
  double t2211;
  double t2213;
  double t2214;
  double t2215;
  double t2217;
  double t2218;
  double t2219;
  double t2374;
  double t2375;
  double t2378;
  double t2380;
  double t2381;
  double t2399;
  double t2400;
  double t2401;
  double t2018;
  double t2035;
  double t2036;
  double t2038;
  double t2039;
  double t2041;
  double t2042;
  double t2045;
  double t2046;
  double t2051;
  double t2084;
  double t2085;
  double t2128;
  double t2129;
  double t2132;
  double t2429;
  double t2431;
  double t2432;
  double t2433;
  double t2256;
  double t2257;
  double t2258;
  double t2259;
  double t2265;
  double t2266;
  double t2268;
  double t2269;
  double t2280;
  double t2281;
  double t2296;
  double t2298;
  double t2299;
  double t2445;
  double t2462;
  double t2463;
  double t2465;
  t1693 = Cos(var1[3]);
  t1769 = Sin(var1[2]);
  t1622 = Cos(var1[2]);
  t1770 = Sin(var1[3]);
  t1780 = Cos(var1[4]);
  t1791 = Sin(var1[4]);
  t1761 = t1622*t1693;
  t1775 = t1769*t1770;
  t1779 = t1761 + t1775;
  t1806 = -1.*t1693*t1769;
  t1807 = t1622*t1770;
  t1822 = t1806 + t1807;
  t1888 = 0.0265*t1693;
  t1894 = -0.0695*t1770;
  t1895 = t1888 + t1894;
  t1897 = t1622*t1895;
  t1898 = -0.0695*t1693;
  t1901 = -0.0265*t1770;
  t1902 = t1898 + t1901;
  t1911 = -1.*t1769*t1902;
  t1797 = -0.2375*t1791;
  t1916 = t1693*t1769;
  t1918 = -1.*t1622*t1770;
  t1919 = t1916 + t1918;
  t1928 = -1.*t1780;
  t1929 = 1. + t1928;
  t1842 = t1780*t1779;
  t1857 = -1.*t1779*t1791;
  t1987 = Cos(var1[5]);
  t1989 = Sin(var1[5]);
  t1992 = Cos(var1[6]);
  t1999 = Sin(var1[6]);
  t1988 = t1622*t1987;
  t1990 = -1.*t1769*t1989;
  t1991 = t1988 + t1990;
  t2019 = -1.*t1987*t1769;
  t2022 = -1.*t1622*t1989;
  t2027 = t2019 + t2022;
  t2056 = -0.0265*t1987;
  t2057 = -0.0695*t1989;
  t2058 = t2056 + t2057;
  t2060 = t1622*t2058;
  t2061 = 0.0695*t1987;
  t2062 = -0.0265*t1989;
  t2065 = t2061 + t2062;
  t2069 = -1.*t1769*t2065;
  t2032 = -0.0265*t1999;
  t2075 = -1.*t1622*t1987;
  t2077 = t1769*t1989;
  t2079 = t2075 + t2077;
  t2086 = -1.*t1992;
  t2088 = 1. + t2086;
  t2037 = -1.*t2027*t1999;
  t2040 = t1992*t2027;
  t2138 = -1.*t1693;
  t2139 = 1. + t2138;
  t2142 = -0.0265*t2139;
  t2143 = t2142 + t1894;
  t2144 = -1.*t1622*t2143;
  t2145 = -0.0695*t2139;
  t2146 = 0.0265*t1770;
  t2147 = t2145 + t2146;
  t2148 = -1.*t1769*t2147;
  t2154 = -1.*t1622*t1693;
  t2155 = -1.*t1769*t1770;
  t2156 = t2154 + t2155;
  t1932 = -0.0265*t1929;
  t1938 = t1932 + t1797;
  t1947 = -0.2375*t1929;
  t1949 = 0.0265*t1791;
  t1950 = t1947 + t1949;
  t2180 = -1.*t1987;
  t2183 = 1. + t2180;
  t2074 = -0.15232*t2027;
  t2082 = -0.025413*t2079;
  t2184 = -0.0695*t2183;
  t2186 = t2184 + t2062;
  t2190 = -1.*t1769*t2186;
  t2192 = -0.0265*t2183;
  t2197 = 0.0695*t1989;
  t2199 = t2192 + t2197;
  t2200 = -1.*t1622*t2199;
  t2089 = -0.2375*t2088;
  t2090 = t2089 + t2032;
  t2093 = t2027*t2090;
  t2094 = -0.0265*t2088;
  t2099 = 0.2375*t1999;
  t2101 = t2094 + t2099;
  t2104 = t2079*t2101;
  t2107 = t1992*t2079;
  t2108 = t2107 + t2037;
  t2115 = -0.025157*t2108;
  t2119 = t2079*t1999;
  t2122 = t2040 + t2119;
  t2123 = -0.312334*t2122;
  t1789 = 0.0265*t1780;
  t1801 = t1789 + t1797;
  t1828 = -0.2375*t1780;
  t1833 = -0.0265*t1791;
  t1834 = t1828 + t1833;
  t1848 = -1.*t1822*t1791;
  t2169 = t1780*t1822;
  t2170 = -1.*t2156*t1791;
  t2172 = t2169 + t2170;
  t2224 = -1.*t1769*t1895;
  t2225 = -1.*t1622*t1902;
  t1851 = t1842 + t1848;
  t1995 = -0.0265*t1992;
  t2000 = -0.2375*t1999;
  t2002 = t1995 + t2000;
  t2031 = 0.2375*t1992;
  t2034 = t2031 + t2032;
  t2270 = -1.*t1769*t2058;
  t2271 = -1.*t1622*t2065;
  t2272 = t1987*t1769;
  t2273 = t1622*t1989;
  t2275 = t2272 + t2273;
  t2260 = -1.*t2079*t1999;
  t2309 = t1769*t2143;
  t2311 = -1.*t1622*t2147;
  t2164 = t1780*t2156;
  t1976 = t1780*t1919;
  t2277 = -0.025413*t2275;
  t2279 = -0.15232*t2079;
  t2328 = -1.*t1622*t2186;
  t2329 = t1769*t2199;
  t2284 = t2079*t2090;
  t2285 = t2275*t2101;
  t2287 = t2275*t1999;
  t2288 = t2107 + t2287;
  t2289 = -0.312334*t2288;
  t2293 = t1992*t2275;
  t2294 = t2293 + t2260;
  t2295 = -0.025157*t2294;
  t2318 = -1.*t1919*t1791;
  t2319 = t2164 + t2318;
  t2151 = -0.15232*t1822;
  t2157 = -0.025367*t2156;
  t2350 = t1769*t1902;
  t2351 = -0.0265*t1693;
  t2352 = 0.0695*t1770;
  t2355 = t2351 + t2352;
  t2356 = t1622*t2355;
  t2161 = t2156*t1938;
  t2162 = t1822*t1950;
  t2165 = t1822*t1791;
  t2166 = t2164 + t2165;
  t2167 = -0.025159*t2166;
  t2174 = -0.312342*t2172;
  t1922 = -0.15232*t1919;
  t1925 = -0.025367*t1779;
  t1926 = t1897 + t1911 + t1922 + t1925;
  t1927 = 0.69051*t1926;
  t1943 = t1779*t1938;
  t1952 = t1919*t1950;
  t1953 = t1919*t1791;
  t1966 = t1842 + t1953;
  t1970 = -0.025159*t1966;
  t1977 = t1976 + t1857;
  t1979 = -0.312342*t1977;
  t1980 = t1897 + t1911 + t1943 + t1952 + t1970 + t1979;
  t1981 = 0.19605*t1980;
  t1982 = t1927 + t1981;
  t2312 = -0.025367*t1919;
  t2313 = -0.15232*t2156;
  t2383 = t1622*t1902;
  t2384 = -1.*t1769*t2355;
  t2316 = t1919*t1938;
  t2317 = t2156*t1950;
  t2320 = -0.312342*t2319;
  t2323 = t2156*t1791;
  t2324 = t1976 + t2323;
  t2325 = -0.025159*t2324;
  t2229 = -0.025367*t1822;
  t2230 = -0.15232*t1779;
  t2231 = t2224 + t2225 + t2229 + t2230;
  t2233 = 0.69051*t2231;
  t2234 = t1822*t1938;
  t2235 = t1779*t1950;
  t2236 = -0.312342*t1851;
  t2237 = t1779*t1791;
  t2238 = t2169 + t2237;
  t2243 = -0.025159*t2238;
  t2246 = t2224 + t2225 + t2234 + t2235 + t2236 + t2243;
  t2247 = 0.19605*t2246;
  t2248 = t2233 + t2247;
  t2339 = t2156*t1801;
  t2341 = t1919*t1834;
  t2342 = -0.025159*t2319;
  t2343 = -1.*t1780*t1919;
  t2345 = t2343 + t2170;
  t2346 = -0.312342*t2345;
  t2347 = t2339 + t2341 + t2342 + t2346;
  t1803 = t1779*t1801;
  t1836 = t1822*t1834;
  t1855 = -0.025159*t1851;
  t1856 = -1.*t1780*t1822;
  t1865 = t1856 + t1857;
  t1869 = -0.312342*t1865;
  t1884 = t1803 + t1836 + t1855 + t1869;
  t2376 = -1.*t1780*t1779;
  t2377 = t2376 + t2318;
  t2211 = t1822*t1801;
  t2213 = t2156*t1834;
  t2214 = -1.*t1780*t2156;
  t2215 = t2214 + t1848;
  t2217 = -0.312342*t2215;
  t2218 = -0.025159*t2172;
  t2219 = t2211 + t2213 + t2217 + t2218;
  t2374 = t1919*t1801;
  t2375 = t1779*t1834;
  t2378 = -0.312342*t2377;
  t2380 = -0.025159*t1977;
  t2381 = t2374 + t2375 + t2378 + t2380;
  t2399 = -0.0265*t1780;
  t2400 = 0.2375*t1791;
  t2401 = t2399 + t2400;
  t2018 = t1991*t2002;
  t2035 = t2027*t2034;
  t2036 = -1.*t1992*t1991;
  t2038 = t2036 + t2037;
  t2039 = -0.025157*t2038;
  t2041 = -1.*t1991*t1999;
  t2042 = t2040 + t2041;
  t2045 = -0.312334*t2042;
  t2046 = t2018 + t2035 + t2039 + t2045;
  t2051 = 0.05554264927529662*var2[6]*t2046;
  t2084 = t2060 + t2069 + t2074 + t2082;
  t2085 = 0.69051*t2084;
  t2128 = t2060 + t2069 + t2093 + t2104 + t2115 + t2123;
  t2129 = 0.19605*t2128;
  t2132 = t2085 + t2129;
  t2429 = -0.0695*t1987;
  t2431 = 0.0265*t1989;
  t2432 = t2429 + t2431;
  t2433 = t1769*t2432;
  t2256 = t2027*t2002;
  t2257 = t2079*t2034;
  t2258 = -0.312334*t2108;
  t2259 = -1.*t1992*t2027;
  t2265 = t2259 + t2260;
  t2266 = -0.025157*t2265;
  t2268 = t2256 + t2257 + t2258 + t2266;
  t2269 = 0.05554264927529662*var2[6]*t2268;
  t2280 = t2270 + t2271 + t2277 + t2279;
  t2281 = 0.69051*t2280;
  t2296 = t2270 + t2271 + t2284 + t2285 + t2289 + t2295;
  t2298 = 0.19605*t2296;
  t2299 = t2281 + t2298;
  t2445 = t1622*t2432;
  t2462 = -0.2375*t1992;
  t2463 = 0.0265*t1999;
  t2465 = t2462 + t2463;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=t2051 + 0.28330859104971495*(1.7566*(0.026461*t1622 + 0.046589*t1769) + 0.69051*(t2144 + t2148 + t2151 + t2157) + 0.19605*(t2144 + t2148 + t2161 + t2162 + t2167 + t2174) + 0.69051*(t2074 + t2082 + t2190 + t2200) + 0.19605*(t2093 + t2104 + t2115 + t2123 + t2190 + t2200))*var2[2] + 0.28330859104971495*t1982*var2[3] + 0.05554264927529662*t1884*var2[4] + 0.28330859104971495*t2132*var2[5];
  p_output1[7]=0;
  p_output1[8]=t2269 + 0.28330859104971495*(1.7566*(0.046589*t1622 - 0.026461*t1769) + 0.69051*(t2309 + t2311 + t2312 + t2313) + 0.19605*(t2309 + t2311 + t2316 + t2317 + t2320 + t2325) + 0.69051*(t2277 + t2279 + t2328 + t2329) + 0.19605*(t2284 + t2285 + t2289 + t2295 + t2328 + t2329))*var2[2] + 0.28330859104971495*t2248*var2[3] + 0.05554264927529662*t2219*var2[4] + 0.28330859104971495*t2299*var2[5];
  p_output1[9]=0.28330859104971495*t1982*var2[2] + 0.28330859104971495*(0.69051*(t2151 + t2157 + t2350 + t2356) + 0.19605*(t2161 + t2162 + t2167 + t2174 + t2350 + t2356))*var2[3] + 0.05554264927529662*t2347*var2[4];
  p_output1[10]=0;
  p_output1[11]=0.28330859104971495*t2248*var2[2] + 0.28330859104971495*(0.69051*(t2312 + t2313 + t2383 + t2384) + 0.19605*(t2316 + t2317 + t2320 + t2325 + t2383 + t2384))*var2[3] + 0.05554264927529662*t2381*var2[4];
  p_output1[12]=0.05554264927529662*t1884*var2[2] + 0.05554264927529662*t2347*var2[3] + 0.05554264927529662*(t2341 - 0.312342*(t2237 + t2343) - 0.025159*t2377 + t1779*t2401)*var2[4];
  p_output1[13]=0;
  p_output1[14]=0.05554264927529662*t2219*var2[2] + 0.05554264927529662*t2381*var2[3] + 0.05554264927529662*(-0.025159*t1865 + t2375 - 0.312342*(t2165 + t2376) + t1822*t2401)*var2[4];
  p_output1[15]=t2051 + 0.28330859104971495*t2132*var2[2] + 0.28330859104971495*(0.69051*(t2060 + t2074 + t2082 + t2433) + 0.19605*(t2060 + t2093 + t2104 + t2115 + t2123 + t2433))*var2[5];
  p_output1[16]=0;
  p_output1[17]=t2269 + 0.28330859104971495*t2299*var2[2] + 0.28330859104971495*(0.69051*(t2270 + t2277 + t2279 + t2445) + 0.19605*(t2270 + t2284 + t2285 + t2289 + t2295 + t2445))*var2[5];
  p_output1[18]=0.05554264927529662*t2046*var2[2] + 0.05554264927529662*t2046*var2[5] + 0.05554264927529662*(t2018 - 0.312334*(t2041 - 1.*t1992*t2275) - 0.025157*(t2036 + t2287) + t2275*t2465)*var2[6];
  p_output1[19]=0;
  p_output1[20]=0.05554264927529662*t2268*var2[2] + 0.05554264927529662*t2268*var2[5] + 0.05554264927529662*(-0.312334*t2038 + t2256 - 0.025157*(t1991*t1999 + t2259) + t1991*t2465)*var2[6];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "dJ_COMPosition.hh"

namespace SymFunction
{

void dJ_COMPosition_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
