/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:44 GMT-05:00
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
  double t1976;
  double t2039;
  double t2089;
  double t134;
  double t2146;
  double t2151;
  double t2158;
  double t2162;
  double t2134;
  double t2136;
  double t2142;
  double t2172;
  double t2174;
  double t2176;
  double t2123;
  double t2223;
  double t2229;
  double t2159;
  double t2164;
  double t2166;
  double t2169;
  double t2179;
  double t2203;
  double t2204;
  double t2205;
  double t2206;
  double t2207;
  double t2208;
  double t2209;
  double t2211;
  double t2213;
  double t2215;
  double t2217;
  double t2256;
  double t2257;
  double t2258;
  double t2265;
  double t2266;
  double t2268;
  double t2272;
  double t2273;
  double t1982;
  double t2057;
  double t2061;
  double t2115;
  double t2129;
  double t2314;
  double t2315;
  double t2316;
  double t2299;
  double t2231;
  double t2233;
  double t2237;
  double t2243;
  double t2246;
  double t2313;
  double t2317;
  double t2319;
  double t2320;
  double t2323;
  double t2324;
  double t2325;
  double t2326;
  t1976 = Cos(var1[5]);
  t2039 = Sin(var1[5]);
  t2089 = Cos(var1[2]);
  t134 = Sin(var1[2]);
  t2146 = Cos(var1[6]);
  t2151 = -1.*t2146;
  t2158 = 1. + t2151;
  t2162 = Sin(var1[6]);
  t2134 = t2089*t1976;
  t2136 = -1.*t134*t2039;
  t2142 = t2134 + t2136;
  t2172 = -1.*t1976*t134;
  t2174 = -1.*t2089*t2039;
  t2176 = t2172 + t2174;
  t2123 = -0.0265*t2039;
  t2223 = -1.*t1976;
  t2229 = 1. + t2223;
  t2159 = -0.2375*t2158;
  t2164 = -0.0265*t2162;
  t2166 = t2159 + t2164;
  t2169 = t2142*t2166;
  t2179 = -0.0265*t2158;
  t2203 = 0.2375*t2162;
  t2204 = t2179 + t2203;
  t2205 = t2176*t2204;
  t2206 = t2146*t2142;
  t2207 = t2176*t2162;
  t2208 = t2206 + t2207;
  t2209 = 0.0325*t2208;
  t2211 = t2146*t2176;
  t2213 = -1.*t2142*t2162;
  t2215 = t2211 + t2213;
  t2217 = -0.0265*t2215;
  t2256 = t1976*t134;
  t2257 = t2089*t2039;
  t2258 = t2256 + t2257;
  t2265 = -0.0265*t2146;
  t2266 = -0.2375*t2162;
  t2268 = t2265 + t2266;
  t2272 = 0.2375*t2146;
  t2273 = t2272 + t2164;
  t1982 = -0.0265*t1976;
  t2057 = -0.0695*t2039;
  t2061 = t1982 + t2057;
  t2115 = 0.0695*t1976;
  t2129 = t2115 + t2123;
  t2314 = -1.*t2089*t1976;
  t2315 = t134*t2039;
  t2316 = t2314 + t2315;
  t2299 = -1.*t2176*t2162;
  t2231 = -0.0695*t2229;
  t2233 = t2231 + t2123;
  t2237 = -0.0265*t2229;
  t2243 = 0.0695*t2039;
  t2246 = t2237 + t2243;
  t2313 = t2176*t2166;
  t2317 = t2316*t2204;
  t2319 = t2146*t2316;
  t2320 = t2319 + t2299;
  t2323 = -0.0265*t2320;
  t2324 = t2316*t2162;
  t2325 = t2211 + t2324;
  t2326 = 0.0325*t2325;
  p_output1[0]=var1[7] + (t2169 + t2205 + t2209 + t2217 + t2089*t2233 - 1.*t134*t2246)*var1[9] + (t134*t2061 + t2089*t2129 + t2169 + t2205 + t2209 + t2217)*var1[12] + (-0.0265*(t2213 - 1.*t2146*t2258) + 0.0325*(t2206 - 1.*t2162*t2258) + t2258*t2268 + t2142*t2273)*var1[13];
  p_output1[1]=var1[8] + (-1.*t134*t2233 - 1.*t2089*t2246 + t2313 + t2317 + t2323 + t2326)*var1[9] + (t2061*t2089 - 1.*t134*t2129 + t2313 + t2317 + t2323 + t2326)*var1[12] + (0.0325*t2215 + t2142*t2268 + t2176*t2273 - 0.0265*(-1.*t2142*t2146 + t2299))*var1[13];
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
    ( !(mrows == 14 && ncols == 1) && 
      !(mrows == 1 && ncols == 14))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "vLeftToe.hh"

namespace SymFunction
{

void vLeftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
