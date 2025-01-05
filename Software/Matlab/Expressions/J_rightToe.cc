/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:38 GMT-05:00
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
  double t1973;
  double t1982;
  double t1999;
  double t2033;
  double t77;
  double t2096;
  double t2150;
  double t2156;
  double t2158;
  double t2175;
  double t2190;
  double t2191;
  double t2193;
  double t2139;
  double t2141;
  double t2146;
  double t2027;
  double t2036;
  double t2041;
  double t2097;
  double t2114;
  double t2118;
  double t2161;
  double t2176;
  double t2184;
  double t2194;
  double t2199;
  double t2203;
  double t2236;
  double t2243;
  double t2249;
  double t2221;
  double t2261;
  double t2288;
  double t2290;
  double t2291;
  double t2272;
  double t2274;
  double t2276;
  double t2280;
  double t2286;
  double t2215;
  double t2303;
  double t2294;
  double t2330;
  double t2334;
  double t2338;
  double t2341;
  double t2344;
  double t2345;
  double t2346;
  double t2216;
  double t2218;
  t1973 = Cos(var1[3]);
  t1982 = -1.*t1973;
  t1999 = 1. + t1982;
  t2033 = Sin(var1[3]);
  t77 = Sin(var1[2]);
  t2096 = Cos(var1[2]);
  t2150 = Cos(var1[4]);
  t2156 = -1.*t2150;
  t2158 = 1. + t2156;
  t2175 = Sin(var1[4]);
  t2190 = t2096*t1973;
  t2191 = t77*t2033;
  t2193 = t2190 + t2191;
  t2139 = -1.*t1973*t77;
  t2141 = t2096*t2033;
  t2146 = t2139 + t2141;
  t2027 = -0.0265*t1999;
  t2036 = -0.0695*t2033;
  t2041 = t2027 + t2036;
  t2097 = -0.0695*t1999;
  t2114 = 0.0265*t2033;
  t2118 = t2097 + t2114;
  t2161 = -0.0265*t2158;
  t2176 = -0.2375*t2175;
  t2184 = t2161 + t2176;
  t2194 = -0.2375*t2158;
  t2199 = 0.0265*t2175;
  t2203 = t2194 + t2199;
  t2236 = -1.*t2096*t1973;
  t2243 = -1.*t77*t2033;
  t2249 = t2236 + t2243;
  t2221 = t2150*t2146;
  t2261 = t2150*t2249;
  t2288 = t1973*t77;
  t2290 = -1.*t2096*t2033;
  t2291 = t2288 + t2290;
  t2272 = 0.0265*t1973;
  t2274 = t2272 + t2036;
  t2276 = -0.0695*t1973;
  t2280 = -0.0265*t2033;
  t2286 = t2276 + t2280;
  t2215 = t2150*t2193;
  t2303 = t2150*t2291;
  t2294 = -1.*t2291*t2175;
  t2330 = -1.*t2193*t2175;
  t2334 = t2303 + t2330;
  t2338 = 0.0265*t2150;
  t2341 = t2338 + t2176;
  t2344 = -0.2375*t2150;
  t2345 = -0.0265*t2175;
  t2346 = t2344 + t2345;
  t2216 = -1.*t2146*t2175;
  t2218 = t2215 + t2216;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=1.;
  p_output1[4]=t2096*t2118 + t2146*t2184 + t2193*t2203 - 0.0115*t2218 - 0.0265*(t2175*t2193 + t2221) - 1.*t2041*t77;
  p_output1[5]=-1.*t2041*t2096 + t2146*t2203 + t2184*t2249 - 0.0115*(t2221 - 1.*t2175*t2249) - 0.0265*(t2146*t2175 + t2261) - 1.*t2118*t77;
  p_output1[6]=t2203*t2249 + t2096*t2286 + t2184*t2291 - 0.0115*(t2261 + t2294) - 0.0265*(t2175*t2249 + t2303) + t2274*t77;
  p_output1[7]=t2184*t2193 + t2096*t2274 + t2203*t2291 - 0.0265*(t2215 + t2175*t2291) - 0.0115*t2334 - 1.*t2286*t77;
  p_output1[8]=-0.0115*(-1.*t2150*t2193 + t2294) - 0.0265*t2334 + t2291*t2341 + t2193*t2346;
  p_output1[9]=-0.0265*t2218 - 0.0115*(-1.*t2146*t2150 + t2330) + t2193*t2341 + t2146*t2346;
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
