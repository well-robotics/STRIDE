/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:35 GMT-05:00
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
  double t1912;
  double t1939;
  double t1958;
  double t1979;
  double t1893;
  double t2013;
  double t2075;
  double t2089;
  double t2091;
  double t2095;
  double t2041;
  double t2058;
  double t2064;
  double t2102;
  double t2103;
  double t2108;
  double t1973;
  double t1982;
  double t1994;
  double t2027;
  double t2033;
  double t2035;
  double t2092;
  double t2096;
  double t2097;
  double t2114;
  double t2118;
  double t2119;
  double t2170;
  double t2175;
  double t2176;
  double t2146;
  double t2099;
  double t2131;
  double t2135;
  double t2139;
  double t2141;
  double t2142;
  double t2149;
  double t2150;
  double t2154;
  double t2194;
  double t2199;
  double t2202;
  double t2208;
  double t2214;
  double t2164;
  double t2179;
  double t2184;
  double t2185;
  double t2186;
  double t2189;
  double t2190;
  double t2191;
  double t2192;
  double t2222;
  double t2223;
  double t2224;
  double t2226;
  double t2231;
  double t2232;
  double t2236;
  double t2243;
  t1912 = Cos(var1[5]);
  t1939 = -1.*t1912;
  t1958 = 1. + t1939;
  t1979 = Sin(var1[5]);
  t1893 = Cos(var1[2]);
  t2013 = Sin(var1[2]);
  t2075 = Cos(var1[6]);
  t2089 = -1.*t2075;
  t2091 = 1. + t2089;
  t2095 = Sin(var1[6]);
  t2041 = t1893*t1912;
  t2058 = -1.*t2013*t1979;
  t2064 = t2041 + t2058;
  t2102 = -1.*t1912*t2013;
  t2103 = -1.*t1893*t1979;
  t2108 = t2102 + t2103;
  t1973 = -0.0695*t1958;
  t1982 = -0.0265*t1979;
  t1994 = t1973 + t1982;
  t2027 = -0.0265*t1958;
  t2033 = 0.0695*t1979;
  t2035 = t2027 + t2033;
  t2092 = -0.2375*t2091;
  t2096 = -0.0265*t2095;
  t2097 = t2092 + t2096;
  t2114 = -0.0265*t2091;
  t2118 = 0.2375*t2095;
  t2119 = t2114 + t2118;
  t2170 = -1.*t1893*t1912;
  t2175 = t2013*t1979;
  t2176 = t2170 + t2175;
  t2146 = t2075*t2108;
  t2099 = t2064*t2097;
  t2131 = t2108*t2119;
  t2135 = t2075*t2064;
  t2139 = t2108*t2095;
  t2141 = t2135 + t2139;
  t2142 = -0.0115*t2141;
  t2149 = -1.*t2064*t2095;
  t2150 = t2146 + t2149;
  t2154 = -0.0265*t2150;
  t2194 = -0.0265*t1912;
  t2199 = -0.0695*t1979;
  t2202 = t2194 + t2199;
  t2208 = 0.0695*t1912;
  t2214 = t2208 + t1982;
  t2164 = t2108*t2097;
  t2179 = t2176*t2119;
  t2184 = t2075*t2176;
  t2185 = -1.*t2108*t2095;
  t2186 = t2184 + t2185;
  t2189 = -0.0265*t2186;
  t2190 = t2176*t2095;
  t2191 = t2146 + t2190;
  t2192 = -0.0115*t2191;
  t2222 = t1912*t2013;
  t2223 = t1893*t1979;
  t2224 = t2222 + t2223;
  t2226 = -0.0265*t2075;
  t2231 = -0.2375*t2095;
  t2232 = t2226 + t2231;
  t2236 = 0.2375*t2075;
  t2243 = t2236 + t2096;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=1.;
  p_output1[4]=t1893*t1994 - 1.*t2013*t2035 + t2099 + t2131 + t2142 + t2154;
  p_output1[5]=-1.*t1994*t2013 - 1.*t1893*t2035 + t2164 + t2179 + t2189 + t2192;
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=t2099 + t2131 + t2142 + t2154 + t2013*t2202 + t1893*t2214;
  p_output1[11]=t2164 + t2179 + t2189 + t2192 + t1893*t2202 - 1.*t2013*t2214;
  p_output1[12]=-0.0265*(t2149 - 1.*t2075*t2224) - 0.0115*(t2135 - 1.*t2095*t2224) + t2224*t2232 + t2064*t2243;
  p_output1[13]=-0.0115*t2150 - 0.0265*(-1.*t2064*t2075 + t2185) + t2064*t2232 + t2108*t2243;
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

#include "J_leftToe.hh"

namespace SymFunction
{

void J_leftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
