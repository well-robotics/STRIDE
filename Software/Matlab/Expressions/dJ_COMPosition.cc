/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:28 GMT-05:00
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
  double t1802;
  double t1826;
  double t1746;
  double t1827;
  double t1844;
  double t1848;
  double t1825;
  double t1832;
  double t1836;
  double t1863;
  double t1864;
  double t1882;
  double t1950;
  double t1951;
  double t1952;
  double t1954;
  double t1955;
  double t1958;
  double t1959;
  double t1968;
  double t1854;
  double t1973;
  double t1975;
  double t1976;
  double t1985;
  double t1986;
  double t1899;
  double t1921;
  double t2044;
  double t2046;
  double t2049;
  double t2056;
  double t2045;
  double t2047;
  double t2048;
  double t2076;
  double t2079;
  double t2086;
  double t2113;
  double t2114;
  double t2115;
  double t2117;
  double t2118;
  double t2119;
  double t2124;
  double t2128;
  double t2089;
  double t2133;
  double t2135;
  double t2137;
  double t2143;
  double t2145;
  double t2094;
  double t2097;
  double t2195;
  double t2196;
  double t2199;
  double t2200;
  double t2201;
  double t2202;
  double t2203;
  double t2205;
  double t2206;
  double t2211;
  double t2212;
  double t2213;
  double t1994;
  double t1999;
  double t2005;
  double t2006;
  double t2007;
  double t2237;
  double t2240;
  double t2131;
  double t2139;
  double t2243;
  double t2247;
  double t2248;
  double t2249;
  double t2254;
  double t2256;
  double t2257;
  double t2146;
  double t2149;
  double t2150;
  double t2154;
  double t2156;
  double t2158;
  double t2161;
  double t2164;
  double t2170;
  double t2175;
  double t2176;
  double t2179;
  double t2184;
  double t1846;
  double t1858;
  double t1885;
  double t1890;
  double t1891;
  double t1905;
  double t2226;
  double t2227;
  double t2231;
  double t2283;
  double t2285;
  double t1908;
  double t2052;
  double t2058;
  double t2064;
  double t2088;
  double t2091;
  double t2327;
  double t2328;
  double t2329;
  double t2330;
  double t2332;
  double t2317;
  double t2366;
  double t2368;
  double t2221;
  double t2033;
  double t2334;
  double t2336;
  double t2385;
  double t2386;
  double t2341;
  double t2342;
  double t2344;
  double t2345;
  double t2346;
  double t2350;
  double t2351;
  double t2352;
  double t2375;
  double t2376;
  double t2208;
  double t2214;
  double t2407;
  double t2408;
  double t2409;
  double t2412;
  double t2413;
  double t2218;
  double t2219;
  double t2222;
  double t2223;
  double t2224;
  double t2232;
  double t1979;
  double t1982;
  double t1983;
  double t1984;
  double t2002;
  double t2009;
  double t2013;
  double t2023;
  double t2027;
  double t2035;
  double t2036;
  double t2037;
  double t2038;
  double t2041;
  double t2369;
  double t2370;
  double t2440;
  double t2441;
  double t2373;
  double t2374;
  double t2377;
  double t2380;
  double t2381;
  double t2382;
  double t2286;
  double t2287;
  double t2288;
  double t2290;
  double t2291;
  double t2292;
  double t2293;
  double t2294;
  double t2295;
  double t2300;
  double t2303;
  double t2304;
  double t2308;
  double t2396;
  double t2398;
  double t2399;
  double t2400;
  double t2402;
  double t2403;
  double t2404;
  double t1861;
  double t1893;
  double t1912;
  double t1913;
  double t1922;
  double t1939;
  double t1941;
  double t2433;
  double t2434;
  double t2268;
  double t2270;
  double t2271;
  double t2272;
  double t2274;
  double t2275;
  double t2276;
  double t2431;
  double t2432;
  double t2435;
  double t2437;
  double t2438;
  double t2456;
  double t2457;
  double t2458;
  double t2075;
  double t2092;
  double t2093;
  double t2095;
  double t2096;
  double t2098;
  double t2099;
  double t2102;
  double t2103;
  double t2108;
  double t2141;
  double t2142;
  double t2185;
  double t2186;
  double t2189;
  double t2486;
  double t2488;
  double t2489;
  double t2490;
  double t2313;
  double t2314;
  double t2315;
  double t2316;
  double t2322;
  double t2323;
  double t2325;
  double t2326;
  double t2337;
  double t2338;
  double t2353;
  double t2355;
  double t2356;
  double t2502;
  double t2519;
  double t2520;
  double t2522;
  t1802 = Cos(var1[3]);
  t1826 = Sin(var1[2]);
  t1746 = Cos(var1[2]);
  t1827 = Sin(var1[3]);
  t1844 = Cos(var1[4]);
  t1848 = Sin(var1[4]);
  t1825 = t1746*t1802;
  t1832 = t1826*t1827;
  t1836 = t1825 + t1832;
  t1863 = -1.*t1802*t1826;
  t1864 = t1746*t1827;
  t1882 = t1863 + t1864;
  t1950 = 0.0265*t1802;
  t1951 = -0.0695*t1827;
  t1952 = t1950 + t1951;
  t1954 = t1746*t1952;
  t1955 = -0.0695*t1802;
  t1958 = -0.0265*t1827;
  t1959 = t1955 + t1958;
  t1968 = -1.*t1826*t1959;
  t1854 = -0.2375*t1848;
  t1973 = t1802*t1826;
  t1975 = -1.*t1746*t1827;
  t1976 = t1973 + t1975;
  t1985 = -1.*t1844;
  t1986 = 1. + t1985;
  t1899 = t1844*t1836;
  t1921 = -1.*t1836*t1848;
  t2044 = Cos(var1[5]);
  t2046 = Sin(var1[5]);
  t2049 = Cos(var1[6]);
  t2056 = Sin(var1[6]);
  t2045 = t1746*t2044;
  t2047 = -1.*t1826*t2046;
  t2048 = t2045 + t2047;
  t2076 = -1.*t2044*t1826;
  t2079 = -1.*t1746*t2046;
  t2086 = t2076 + t2079;
  t2113 = -0.0265*t2044;
  t2114 = -0.0695*t2046;
  t2115 = t2113 + t2114;
  t2117 = t1746*t2115;
  t2118 = 0.0695*t2044;
  t2119 = -0.0265*t2046;
  t2124 = t2118 + t2119;
  t2128 = -1.*t1826*t2124;
  t2089 = -0.0265*t2056;
  t2133 = -1.*t1746*t2044;
  t2135 = t1826*t2046;
  t2137 = t2133 + t2135;
  t2143 = -1.*t2049;
  t2145 = 1. + t2143;
  t2094 = -1.*t2086*t2056;
  t2097 = t2049*t2086;
  t2195 = -1.*t1802;
  t2196 = 1. + t2195;
  t2199 = -0.0265*t2196;
  t2200 = t2199 + t1951;
  t2201 = -1.*t1746*t2200;
  t2202 = -0.0695*t2196;
  t2203 = 0.0265*t1827;
  t2205 = t2202 + t2203;
  t2206 = -1.*t1826*t2205;
  t2211 = -1.*t1746*t1802;
  t2212 = -1.*t1826*t1827;
  t2213 = t2211 + t2212;
  t1994 = -0.0265*t1986;
  t1999 = t1994 + t1854;
  t2005 = -0.2375*t1986;
  t2006 = 0.0265*t1848;
  t2007 = t2005 + t2006;
  t2237 = -1.*t2044;
  t2240 = 1. + t2237;
  t2131 = -0.15232*t2086;
  t2139 = -0.025413*t2137;
  t2243 = -0.0695*t2240;
  t2247 = t2243 + t2119;
  t2248 = -1.*t1826*t2247;
  t2249 = -0.0265*t2240;
  t2254 = 0.0695*t2046;
  t2256 = t2249 + t2254;
  t2257 = -1.*t1746*t2256;
  t2146 = -0.2375*t2145;
  t2149 = t2146 + t2089;
  t2150 = t2086*t2149;
  t2154 = -0.0265*t2145;
  t2156 = 0.2375*t2056;
  t2158 = t2154 + t2156;
  t2161 = t2137*t2158;
  t2164 = t2049*t2137;
  t2170 = t2164 + t2094;
  t2175 = -0.025226*t2170;
  t2176 = t2137*t2056;
  t2179 = t2097 + t2176;
  t2184 = -0.314506*t2179;
  t1846 = 0.0265*t1844;
  t1858 = t1846 + t1854;
  t1885 = -0.2375*t1844;
  t1890 = -0.0265*t1848;
  t1891 = t1885 + t1890;
  t1905 = -1.*t1882*t1848;
  t2226 = t1844*t1882;
  t2227 = -1.*t2213*t1848;
  t2231 = t2226 + t2227;
  t2283 = -1.*t1826*t1952;
  t2285 = -1.*t1746*t1959;
  t1908 = t1899 + t1905;
  t2052 = -0.0265*t2049;
  t2058 = -0.2375*t2056;
  t2064 = t2052 + t2058;
  t2088 = 0.2375*t2049;
  t2091 = t2088 + t2089;
  t2327 = -1.*t1826*t2115;
  t2328 = -1.*t1746*t2124;
  t2329 = t2044*t1826;
  t2330 = t1746*t2046;
  t2332 = t2329 + t2330;
  t2317 = -1.*t2137*t2056;
  t2366 = t1826*t2200;
  t2368 = -1.*t1746*t2205;
  t2221 = t1844*t2213;
  t2033 = t1844*t1976;
  t2334 = -0.025413*t2332;
  t2336 = -0.15232*t2137;
  t2385 = -1.*t1746*t2247;
  t2386 = t1826*t2256;
  t2341 = t2137*t2149;
  t2342 = t2332*t2158;
  t2344 = t2332*t2056;
  t2345 = t2164 + t2344;
  t2346 = -0.314506*t2345;
  t2350 = t2049*t2332;
  t2351 = t2350 + t2317;
  t2352 = -0.025226*t2351;
  t2375 = -1.*t1976*t1848;
  t2376 = t2221 + t2375;
  t2208 = -0.15232*t1882;
  t2214 = -0.025367*t2213;
  t2407 = t1826*t1959;
  t2408 = -0.0265*t1802;
  t2409 = 0.0695*t1827;
  t2412 = t2408 + t2409;
  t2413 = t1746*t2412;
  t2218 = t2213*t1999;
  t2219 = t1882*t2007;
  t2222 = t1882*t1848;
  t2223 = t2221 + t2222;
  t2224 = -0.025229*t2223;
  t2232 = -0.314514*t2231;
  t1979 = -0.15232*t1976;
  t1982 = -0.025367*t1836;
  t1983 = t1954 + t1968 + t1979 + t1982;
  t1984 = 0.69051*t1983;
  t2002 = t1836*t1999;
  t2009 = t1976*t2007;
  t2013 = t1976*t1848;
  t2023 = t1899 + t2013;
  t2027 = -0.025229*t2023;
  t2035 = t2033 + t1921;
  t2036 = -0.314514*t2035;
  t2037 = t1954 + t1968 + t2002 + t2009 + t2027 + t2036;
  t2038 = 0.19964*t2037;
  t2041 = t1984 + t2038;
  t2369 = -0.025367*t1976;
  t2370 = -0.15232*t2213;
  t2440 = t1746*t1959;
  t2441 = -1.*t1826*t2412;
  t2373 = t1976*t1999;
  t2374 = t2213*t2007;
  t2377 = -0.314514*t2376;
  t2380 = t2213*t1848;
  t2381 = t2033 + t2380;
  t2382 = -0.025229*t2381;
  t2286 = -0.025367*t1882;
  t2287 = -0.15232*t1836;
  t2288 = t2283 + t2285 + t2286 + t2287;
  t2290 = 0.69051*t2288;
  t2291 = t1882*t1999;
  t2292 = t1836*t2007;
  t2293 = -0.314514*t1908;
  t2294 = t1836*t1848;
  t2295 = t2226 + t2294;
  t2300 = -0.025229*t2295;
  t2303 = t2283 + t2285 + t2291 + t2292 + t2293 + t2300;
  t2304 = 0.19964*t2303;
  t2308 = t2290 + t2304;
  t2396 = t2213*t1858;
  t2398 = t1976*t1891;
  t2399 = -0.025229*t2376;
  t2400 = -1.*t1844*t1976;
  t2402 = t2400 + t2227;
  t2403 = -0.314514*t2402;
  t2404 = t2396 + t2398 + t2399 + t2403;
  t1861 = t1836*t1858;
  t1893 = t1882*t1891;
  t1912 = -0.025229*t1908;
  t1913 = -1.*t1844*t1882;
  t1922 = t1913 + t1921;
  t1939 = -0.314514*t1922;
  t1941 = t1861 + t1893 + t1912 + t1939;
  t2433 = -1.*t1844*t1836;
  t2434 = t2433 + t2375;
  t2268 = t1882*t1858;
  t2270 = t2213*t1891;
  t2271 = -1.*t1844*t2213;
  t2272 = t2271 + t1905;
  t2274 = -0.314514*t2272;
  t2275 = -0.025229*t2231;
  t2276 = t2268 + t2270 + t2274 + t2275;
  t2431 = t1976*t1858;
  t2432 = t1836*t1891;
  t2435 = -0.314514*t2434;
  t2437 = -0.025229*t2035;
  t2438 = t2431 + t2432 + t2435 + t2437;
  t2456 = -0.0265*t1844;
  t2457 = 0.2375*t1848;
  t2458 = t2456 + t2457;
  t2075 = t2048*t2064;
  t2092 = t2086*t2091;
  t2093 = -1.*t2049*t2048;
  t2095 = t2093 + t2094;
  t2096 = -0.025226*t2095;
  t2098 = -1.*t2048*t2056;
  t2099 = t2097 + t2098;
  t2102 = -0.314506*t2099;
  t2103 = t2075 + t2092 + t2096 + t2102;
  t2108 = 0.05982798405705895*var2[6]*t2103;
  t2141 = t2117 + t2128 + t2131 + t2139;
  t2142 = 0.69051*t2141;
  t2185 = t2117 + t2128 + t2150 + t2161 + t2175 + t2184;
  t2186 = 0.19964*t2185;
  t2189 = t2142 + t2186;
  t2486 = -0.0695*t2044;
  t2488 = 0.0265*t2046;
  t2489 = t2486 + t2488;
  t2490 = t1826*t2489;
  t2313 = t2086*t2064;
  t2314 = t2137*t2091;
  t2315 = -0.314506*t2170;
  t2316 = -1.*t2049*t2086;
  t2322 = t2316 + t2317;
  t2323 = -0.025226*t2322;
  t2325 = t2313 + t2314 + t2315 + t2323;
  t2326 = 0.05982798405705895*var2[6]*t2325;
  t2337 = t2327 + t2328 + t2334 + t2336;
  t2338 = 0.69051*t2337;
  t2353 = t2327 + t2328 + t2341 + t2342 + t2346 + t2352;
  t2355 = 0.19964*t2353;
  t2356 = t2338 + t2355;
  t2502 = t1746*t2489;
  t2519 = -0.2375*t2049;
  t2520 = 0.0265*t2056;
  t2522 = t2519 + t2520;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=t2108 + 0.2996793431028799*(1.5566*(0.026461*t1746 + 0.046589*t1826) + 0.69051*(t2201 + t2206 + t2208 + t2214) + 0.19964*(t2201 + t2206 + t2218 + t2219 + t2224 + t2232) + 0.69051*(t2131 + t2139 + t2248 + t2257) + 0.19964*(t2150 + t2161 + t2175 + t2184 + t2248 + t2257))*var2[2] + 0.2996793431028799*t2041*var2[3] + 0.05982798405705895*t1941*var2[4] + 0.2996793431028799*t2189*var2[5];
  p_output1[7]=0;
  p_output1[8]=t2326 + 0.2996793431028799*(1.5566*(0.046589*t1746 - 0.026461*t1826) + 0.69051*(t2366 + t2368 + t2369 + t2370) + 0.19964*(t2366 + t2368 + t2373 + t2374 + t2377 + t2382) + 0.69051*(t2334 + t2336 + t2385 + t2386) + 0.19964*(t2341 + t2342 + t2346 + t2352 + t2385 + t2386))*var2[2] + 0.2996793431028799*t2308*var2[3] + 0.05982798405705895*t2276*var2[4] + 0.2996793431028799*t2356*var2[5];
  p_output1[9]=0.2996793431028799*t2041*var2[2] + 0.2996793431028799*(0.69051*(t2208 + t2214 + t2407 + t2413) + 0.19964*(t2218 + t2219 + t2224 + t2232 + t2407 + t2413))*var2[3] + 0.05982798405705895*t2404*var2[4];
  p_output1[10]=0;
  p_output1[11]=0.2996793431028799*t2308*var2[2] + 0.2996793431028799*(0.69051*(t2369 + t2370 + t2440 + t2441) + 0.19964*(t2373 + t2374 + t2377 + t2382 + t2440 + t2441))*var2[3] + 0.05982798405705895*t2438*var2[4];
  p_output1[12]=0.05982798405705895*t1941*var2[2] + 0.05982798405705895*t2404*var2[3] + 0.05982798405705895*(t2398 - 0.314514*(t2294 + t2400) - 0.025229*t2434 + t1836*t2458)*var2[4];
  p_output1[13]=0;
  p_output1[14]=0.05982798405705895*t2276*var2[2] + 0.05982798405705895*t2438*var2[3] + 0.05982798405705895*(-0.025229*t1922 + t2432 - 0.314514*(t2222 + t2433) + t1882*t2458)*var2[4];
  p_output1[15]=t2108 + 0.2996793431028799*t2189*var2[2] + 0.2996793431028799*(0.69051*(t2117 + t2131 + t2139 + t2490) + 0.19964*(t2117 + t2150 + t2161 + t2175 + t2184 + t2490))*var2[5];
  p_output1[16]=0;
  p_output1[17]=t2326 + 0.2996793431028799*t2356*var2[2] + 0.2996793431028799*(0.69051*(t2327 + t2334 + t2336 + t2502) + 0.19964*(t2327 + t2341 + t2342 + t2346 + t2352 + t2502))*var2[5];
  p_output1[18]=0.05982798405705895*t2103*var2[2] + 0.05982798405705895*t2103*var2[5] + 0.05982798405705895*(t2075 - 0.314506*(t2098 - 1.*t2049*t2332) - 0.025226*(t2093 + t2344) + t2332*t2522)*var2[6];
  p_output1[19]=0;
  p_output1[20]=0.05982798405705895*t2325*var2[2] + 0.05982798405705895*t2325*var2[5] + 0.05982798405705895*(-0.314506*t2095 + t2313 - 0.025226*(t2048*t2056 + t2316) + t2048*t2522)*var2[6];
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
