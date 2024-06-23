/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:38 GMT-05:00
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
  double t1788;
  double t1847;
  double t1906;
  double t1867;
  double t1995;
  double t2012;
  double t2016;
  double t2020;
  double t2032;
  double t2035;
  double t2036;
  double t2040;
  double t2041;
  double t2042;
  double t2083;
  double t2084;
  double t2085;
  double t2043;
  double t2059;
  double t2063;
  double t1996;
  double t2013;
  double t2014;
  double t2086;
  double t2087;
  double t2096;
  double t2123;
  double t2126;
  double t2135;
  double t2137;
  double t2141;
  double t2142;
  double t2143;
  double t2150;
  double t2151;
  double t2152;
  double t2153;
  double t2154;
  double t2159;
  double t2160;
  double t2169;
  double t2155;
  double t2156;
  double t2157;
  double t2136;
  double t2138;
  double t2139;
  double t2170;
  double t2173;
  double t2174;
  double t1938;
  double t1939;
  double t1958;
  double t1970;
  double t1990;
  double t2079;
  double t2097;
  double t2100;
  double t2033;
  double t2117;
  double t2118;
  double t2119;
  double t2129;
  double t2130;
  double t2131;
  double t2132;
  double t2133;
  double t2158;
  double t2181;
  double t2182;
  double t2144;
  double t2194;
  double t2196;
  double t2198;
  double t2233;
  double t2234;
  double t2235;
  double t2236;
  double t2237;
  double t2238;
  double t2239;
  double t2240;
  double t1851;
  double t1923;
  double t1935;
  double t1994;
  double t2211;
  double t2212;
  double t2108;
  double t2253;
  double t2254;
  double t2255;
  double t2207;
  double t2208;
  double t2209;
  double t2250;
  double t2251;
  double t2252;
  double t2256;
  double t2257;
  double t2203;
  double t2204;
  double t2205;
  double t2206;
  double t2269;
  double t2270;
  double t2242;
  double t2243;
  double t2244;
  double t2245;
  double t2288;
  double t2289;
  double t2290;
  double t2291;
  double t2292;
  double t2293;
  double t2294;
  double t2247;
  double t2248;
  double t2249;
  double t2261;
  double t2263;
  double t2264;
  double t2265;
  double t2305;
  double t2306;
  double t2307;
  double t2271;
  double t2273;
  double t2274;
  double t2275;
  double t2276;
  double t2277;
  double t2278;
  double t2325;
  double t2326;
  double t2327;
  double t2328;
  double t2329;
  double t2330;
  double t2331;
  double t2332;
  double t2125;
  double t2127;
  double t2128;
  double t2134;
  double t2226;
  double t2227;
  double t2190;
  double t2345;
  double t2346;
  double t2347;
  double t2220;
  double t2223;
  double t2224;
  double t2342;
  double t2343;
  double t2344;
  double t2348;
  double t2349;
  double t2214;
  double t2215;
  double t2216;
  double t2217;
  double t2361;
  double t2362;
  double t2334;
  double t2335;
  double t2336;
  double t2337;
  double t2380;
  double t2381;
  double t2382;
  double t2383;
  double t2384;
  double t2385;
  double t2386;
  double t2339;
  double t2340;
  double t2341;
  double t2353;
  double t2355;
  double t2356;
  double t2357;
  double t2397;
  double t2398;
  double t2399;
  double t2363;
  double t2365;
  double t2366;
  double t2367;
  double t2368;
  double t2369;
  double t2370;
  t1788 = Cos(var1[2]);
  t1847 = Cos(var1[3]);
  t1906 = Sin(var1[3]);
  t1867 = Sin(var1[2]);
  t1995 = Cos(var1[4]);
  t2012 = Sin(var1[4]);
  t2016 = t1847*t1995;
  t2020 = -1.*t1906*t2012;
  t2032 = t2016 + t2020;
  t2035 = -1.*t1995;
  t2036 = 1. + t2035;
  t2040 = 0.5*t2036;
  t2041 = 0.671885*t1995;
  t2042 = t2040 + t2041;
  t2083 = -1.*t1995*t1906;
  t2084 = -1.*t1847*t2012;
  t2085 = t2083 + t2084;
  t2043 = t2042*t1906;
  t2059 = 0.171885*t1847*t2012;
  t2063 = t2043 + t2059;
  t1996 = t1995*t1906;
  t2013 = t1847*t2012;
  t2014 = t1996 + t2013;
  t2086 = t1847*t2042;
  t2087 = -0.171885*t1906*t2012;
  t2096 = t2086 + t2087;
  t2123 = Cos(var1[5]);
  t2126 = Sin(var1[5]);
  t2135 = Cos(var1[6]);
  t2137 = Sin(var1[6]);
  t2141 = t2123*t2135;
  t2142 = -1.*t2126*t2137;
  t2143 = t2141 + t2142;
  t2150 = -1.*t2135;
  t2151 = 1. + t2150;
  t2152 = 0.5*t2151;
  t2153 = 0.671885*t2135;
  t2154 = t2152 + t2153;
  t2159 = -1.*t2135*t2126;
  t2160 = -1.*t2123*t2137;
  t2169 = t2159 + t2160;
  t2155 = t2154*t2126;
  t2156 = 0.171885*t2123*t2137;
  t2157 = t2155 + t2156;
  t2136 = t2135*t2126;
  t2138 = t2123*t2137;
  t2139 = t2136 + t2138;
  t2170 = t2123*t2154;
  t2173 = -0.171885*t2126*t2137;
  t2174 = t2170 + t2173;
  t1938 = Power(t1847,2);
  t1939 = 0.1494*t1938;
  t1958 = Power(t1906,2);
  t1970 = 0.1494*t1958;
  t1990 = t1939 + t1970;
  t2079 = -1.*t2063*t2032;
  t2097 = -1.*t2085*t2096;
  t2100 = t2079 + t2097;
  t2033 = -1.*t1867*t2032;
  t2117 = t2063*t2014;
  t2118 = t2032*t2096;
  t2119 = t2117 + t2118;
  t2129 = Power(t2123,2);
  t2130 = 0.1494*t2129;
  t2131 = Power(t2126,2);
  t2132 = 0.1494*t2131;
  t2133 = t2130 + t2132;
  t2158 = -1.*t2157*t2143;
  t2181 = -1.*t2169*t2174;
  t2182 = t2158 + t2181;
  t2144 = -1.*t1867*t2143;
  t2194 = t2157*t2139;
  t2196 = t2143*t2174;
  t2198 = t2194 + t2196;
  t2233 = -1.*t2042*t1906;
  t2234 = -0.171885*t1847*t2012;
  t2235 = t2233 + t2234;
  t2236 = t2235*t2032;
  t2237 = t2063*t2032;
  t2238 = t2085*t2096;
  t2239 = t2014*t2096;
  t2240 = t2236 + t2237 + t2238 + t2239;
  t1851 = -1.*t1788*t1847;
  t1923 = t1867*t1906;
  t1935 = t1851 + t1923;
  t1994 = -3.4261*t1935*t1990;
  t2211 = t1788*t2085;
  t2212 = t2211 + t2033;
  t2108 = -1.*t1867*t2085;
  t2253 = -1.*t1847*t1995;
  t2254 = t1906*t2012;
  t2255 = t2253 + t2254;
  t2207 = -1.*t1867*t2014;
  t2208 = t1788*t2032;
  t2209 = t2207 + t2208;
  t2250 = -1.*t2085*t2235;
  t2251 = -1.*t2085*t2063;
  t2252 = -1.*t2032*t2096;
  t2256 = -1.*t2096*t2255;
  t2257 = t2250 + t2251 + t2252 + t2256;
  t2203 = -1.*t1847*t1867;
  t2204 = -1.*t1788*t1906;
  t2205 = t2203 + t2204;
  t2206 = -3.4261*t2205*t1990;
  t2269 = t1867*t2085;
  t2270 = t2269 + t2208;
  t2242 = t2042*t1995;
  t2243 = Power(t2012,2);
  t2244 = 0.171885*t2243;
  t2245 = t2242 + t2244;
  t2288 = -0.171885*t1995*t1906;
  t2289 = t2288 + t2234;
  t2290 = t2289*t2032;
  t2291 = 0.171885*t1847*t1995;
  t2292 = t2291 + t2087;
  t2293 = t2014*t2292;
  t2294 = t2290 + t2237 + t2238 + t2293;
  t2247 = t2042*t2012;
  t2248 = -0.171885*t1995*t2012;
  t2249 = t2247 + t2248;
  t2261 = -0.85216*t2212*t2100;
  t2263 = t1788*t2255;
  t2264 = t2108 + t2263;
  t2265 = -0.85216*t2119*t2264;
  t2305 = -1.*t2085*t2289;
  t2306 = -1.*t2032*t2292;
  t2307 = t2305 + t2251 + t2306 + t2256;
  t2271 = -0.85216*t2270*t2100;
  t2273 = t1867*t2255;
  t2274 = t2211 + t2273;
  t2275 = -0.85216*t2119*t2274;
  t2276 = t1788*t2014;
  t2277 = t1867*t2032;
  t2278 = t2276 + t2277;
  t2325 = -1.*t2154*t2126;
  t2326 = -0.171885*t2123*t2137;
  t2327 = t2325 + t2326;
  t2328 = t2327*t2143;
  t2329 = t2157*t2143;
  t2330 = t2169*t2174;
  t2331 = t2139*t2174;
  t2332 = t2328 + t2329 + t2330 + t2331;
  t2125 = -1.*t1788*t2123;
  t2127 = t1867*t2126;
  t2128 = t2125 + t2127;
  t2134 = -3.4261*t2128*t2133;
  t2226 = t1788*t2169;
  t2227 = t2226 + t2144;
  t2190 = -1.*t1867*t2169;
  t2345 = -1.*t2123*t2135;
  t2346 = t2126*t2137;
  t2347 = t2345 + t2346;
  t2220 = -1.*t1867*t2139;
  t2223 = t1788*t2143;
  t2224 = t2220 + t2223;
  t2342 = -1.*t2169*t2327;
  t2343 = -1.*t2169*t2157;
  t2344 = -1.*t2143*t2174;
  t2348 = -1.*t2174*t2347;
  t2349 = t2342 + t2343 + t2344 + t2348;
  t2214 = -1.*t2123*t1867;
  t2215 = -1.*t1788*t2126;
  t2216 = t2214 + t2215;
  t2217 = -3.4261*t2216*t2133;
  t2361 = t1867*t2169;
  t2362 = t2361 + t2223;
  t2334 = t2154*t2135;
  t2335 = Power(t2137,2);
  t2336 = 0.171885*t2335;
  t2337 = t2334 + t2336;
  t2380 = -0.171885*t2135*t2126;
  t2381 = t2380 + t2326;
  t2382 = t2381*t2143;
  t2383 = 0.171885*t2123*t2135;
  t2384 = t2383 + t2173;
  t2385 = t2139*t2384;
  t2386 = t2382 + t2329 + t2330 + t2385;
  t2339 = t2154*t2137;
  t2340 = -0.171885*t2135*t2137;
  t2341 = t2339 + t2340;
  t2353 = -0.85216*t2227*t2182;
  t2355 = t1788*t2347;
  t2356 = t2190 + t2355;
  t2357 = -0.85216*t2198*t2356;
  t2397 = -1.*t2169*t2381;
  t2398 = -1.*t2143*t2384;
  t2399 = t2397 + t2343 + t2398 + t2348;
  t2363 = -0.85216*t2362*t2182;
  t2365 = t1867*t2347;
  t2366 = t2226 + t2365;
  t2367 = -0.85216*t2198*t2366;
  t2368 = t1788*t2139;
  t2369 = t1867*t2143;
  t2370 = t2368 + t2369;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(-0.5*(3.70591*t1867 + t2206 - 0.85216*t2100*t2209 - 0.85216*t2119*t2212 + t2217 - 0.85216*t2182*t2224 - 0.85216*t2198*t2227)*var2[0] - 0.5*(3.70591*t1788 + t1994 - 0.85216*(-1.*t1788*t2014 + t2033)*t2100 - 0.85216*(-1.*t1788*t2032 + t2108)*t2119 + t2134 - 0.85216*(-1.*t1788*t2139 + t2144)*t2182 - 0.85216*(-1.*t1788*t2143 + t2190)*t2198)*var2[1])*var2[2];
  p_output1[3]=var2[2]*(-0.5*(t2206 - 0.85216*t2240*t2270 + t2271 + t2275 - 0.85216*t2257*t2278)*var2[0] - 0.5*(t1994 - 0.85216*t2212*t2240 - 0.85216*t2209*t2257 + t2261 + t2265)*var2[1] - 0.5*(-1.70432*t2119*t2240 - 1.70432*t2100*t2257)*var2[2] - 0.5*(-0.85216*t2240*t2245 - 0.85216*t2249*t2257)*var2[3] + 0.0732367608*t2240*var2[4]);
  p_output1[4]=var2[2]*(-0.5*(t2271 + t2275 - 0.85216*t2270*t2294 - 0.85216*t2278*t2307)*var2[0] - 0.5*(t2261 + t2265 - 0.85216*t2212*t2294 - 0.85216*t2209*t2307)*var2[1] - 0.5*(-1.70432*t2119*t2294 - 1.70432*t2100*t2307)*var2[2] - 0.5*(-0.85216*(0.171885*t1995*t2012 - 1.*t2012*t2042)*t2119 - 0.85216*t2100*(-0.171885*Power(t1995,2) + t2242) - 0.85216*t2245*t2294 - 0.85216*t2249*t2307)*var2[3] + 0.0732367608*t2294*var2[4]);
  p_output1[5]=var2[2]*(-0.5*(t2217 - 0.85216*t2332*t2362 + t2363 + t2367 - 0.85216*t2349*t2370)*var2[0] - 0.5*(t2134 - 0.85216*t2227*t2332 - 0.85216*t2224*t2349 + t2353 + t2357)*var2[1] - 0.5*(-1.70432*t2198*t2332 - 1.70432*t2182*t2349)*var2[2] - 0.5*(-0.85216*t2332*t2337 - 0.85216*t2341*t2349)*var2[5] + 0.0732367608*t2332*var2[6]);
  p_output1[6]=var2[2]*(-0.5*(t2363 + t2367 - 0.85216*t2362*t2386 - 0.85216*t2370*t2399)*var2[0] - 0.5*(t2353 + t2357 - 0.85216*t2227*t2386 - 0.85216*t2224*t2399)*var2[1] - 0.5*(-1.70432*t2198*t2386 - 1.70432*t2182*t2399)*var2[2] - 0.5*(-0.85216*(0.171885*t2135*t2137 - 1.*t2137*t2154)*t2198 - 0.85216*t2182*(-0.171885*Power(t2135,2) + t2334) - 0.85216*t2337*t2386 - 0.85216*t2341*t2399)*var2[5] + 0.0732367608*t2386*var2[6]);
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

#include "Ce3_vec3_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
