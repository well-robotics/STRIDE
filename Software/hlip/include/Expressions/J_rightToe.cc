/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:43 GMT-05:00
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
  double t1925;
  double t1938;
  double t1970;
  double t1979;
  double t1916;
  double t2040;
  double t2099;
  double t2101;
  double t2104;
  double t2119;
  double t2134;
  double t2136;
  double t2137;
  double t2084;
  double t2089;
  double t2093;
  double t1976;
  double t1982;
  double t2035;
  double t2057;
  double t2061;
  double t2062;
  double t2115;
  double t2123;
  double t2129;
  double t2142;
  double t2146;
  double t2151;
  double t2184;
  double t2192;
  double t2197;
  double t2165;
  double t2205;
  double t2233;
  double t2234;
  double t2235;
  double t2217;
  double t2218;
  double t2223;
  double t2229;
  double t2230;
  double t2159;
  double t2247;
  double t2238;
  double t2277;
  double t2279;
  double t2284;
  double t2285;
  double t2288;
  double t2289;
  double t2293;
  double t2161;
  double t2162;
  t1925 = Cos(var1[3]);
  t1938 = -1.*t1925;
  t1970 = 1. + t1938;
  t1979 = Sin(var1[3]);
  t1916 = Sin(var1[2]);
  t2040 = Cos(var1[2]);
  t2099 = Cos(var1[4]);
  t2101 = -1.*t2099;
  t2104 = 1. + t2101;
  t2119 = Sin(var1[4]);
  t2134 = t2040*t1925;
  t2136 = t1916*t1979;
  t2137 = t2134 + t2136;
  t2084 = -1.*t1925*t1916;
  t2089 = t2040*t1979;
  t2093 = t2084 + t2089;
  t1976 = -0.0265*t1970;
  t1982 = -0.0695*t1979;
  t2035 = t1976 + t1982;
  t2057 = -0.0695*t1970;
  t2061 = 0.0265*t1979;
  t2062 = t2057 + t2061;
  t2115 = -0.0265*t2104;
  t2123 = -0.2375*t2119;
  t2129 = t2115 + t2123;
  t2142 = -0.2375*t2104;
  t2146 = 0.0265*t2119;
  t2151 = t2142 + t2146;
  t2184 = -1.*t2040*t1925;
  t2192 = -1.*t1916*t1979;
  t2197 = t2184 + t2192;
  t2165 = t2099*t2093;
  t2205 = t2099*t2197;
  t2233 = t1925*t1916;
  t2234 = -1.*t2040*t1979;
  t2235 = t2233 + t2234;
  t2217 = 0.0265*t1925;
  t2218 = t2217 + t1982;
  t2223 = -0.0695*t1925;
  t2229 = -0.0265*t1979;
  t2230 = t2223 + t2229;
  t2159 = t2099*t2137;
  t2247 = t2099*t2235;
  t2238 = -1.*t2235*t2119;
  t2277 = -1.*t2137*t2119;
  t2279 = t2247 + t2277;
  t2284 = 0.0265*t2099;
  t2285 = t2284 + t2123;
  t2288 = -0.2375*t2099;
  t2289 = -0.0265*t2119;
  t2293 = t2288 + t2289;
  t2161 = -1.*t2093*t2119;
  t2162 = t2159 + t2161;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=1.;
  p_output1[4]=-1.*t1916*t2035 + t2040*t2062 + t2093*t2129 + t2137*t2151 + 0.0325*t2162 - 0.0265*(t2119*t2137 + t2165);
  p_output1[5]=-1.*t2035*t2040 - 1.*t1916*t2062 + t2093*t2151 + t2129*t2197 + 0.0325*(t2165 - 1.*t2119*t2197) - 0.0265*(t2093*t2119 + t2205);
  p_output1[6]=t2151*t2197 + t1916*t2218 + t2040*t2230 + t2129*t2235 + 0.0325*(t2205 + t2238) - 0.0265*(t2119*t2197 + t2247);
  p_output1[7]=t2129*t2137 + t2040*t2218 - 1.*t1916*t2230 + t2151*t2235 - 0.0265*(t2159 + t2119*t2235) + 0.0325*t2279;
  p_output1[8]=0.0325*(-1.*t2099*t2137 + t2238) - 0.0265*t2279 + t2235*t2285 + t2137*t2293;
  p_output1[9]=-0.0265*t2162 + 0.0325*(-1.*t2093*t2099 + t2277) + t2137*t2285 + t2093*t2293;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "J_rightToe.hh"

namespace SymFunction
{

void J_rightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
