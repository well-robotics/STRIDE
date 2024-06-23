/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:42 GMT-05:00
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
  double t1855;
  double t1869;
  double t1901;
  double t1922;
  double t1836;
  double t1953;
  double t2018;
  double t2032;
  double t2034;
  double t2038;
  double t1982;
  double t2000;
  double t2002;
  double t2045;
  double t2046;
  double t2051;
  double t1916;
  double t1925;
  double t1932;
  double t1970;
  double t1976;
  double t1977;
  double t2035;
  double t2039;
  double t2040;
  double t2057;
  double t2061;
  double t2062;
  double t2108;
  double t2115;
  double t2119;
  double t2089;
  double t2042;
  double t2074;
  double t2077;
  double t2082;
  double t2084;
  double t2085;
  double t2090;
  double t2093;
  double t2094;
  double t2137;
  double t2142;
  double t2145;
  double t2151;
  double t2157;
  double t2107;
  double t2122;
  double t2123;
  double t2128;
  double t2129;
  double t2132;
  double t2133;
  double t2134;
  double t2135;
  double t2165;
  double t2166;
  double t2167;
  double t2169;
  double t2172;
  double t2174;
  double t2179;
  double t2184;
  t1855 = Cos(var1[5]);
  t1869 = -1.*t1855;
  t1901 = 1. + t1869;
  t1922 = Sin(var1[5]);
  t1836 = Cos(var1[2]);
  t1953 = Sin(var1[2]);
  t2018 = Cos(var1[6]);
  t2032 = -1.*t2018;
  t2034 = 1. + t2032;
  t2038 = Sin(var1[6]);
  t1982 = t1836*t1855;
  t2000 = -1.*t1953*t1922;
  t2002 = t1982 + t2000;
  t2045 = -1.*t1855*t1953;
  t2046 = -1.*t1836*t1922;
  t2051 = t2045 + t2046;
  t1916 = -0.0695*t1901;
  t1925 = -0.0265*t1922;
  t1932 = t1916 + t1925;
  t1970 = -0.0265*t1901;
  t1976 = 0.0695*t1922;
  t1977 = t1970 + t1976;
  t2035 = -0.2375*t2034;
  t2039 = -0.0265*t2038;
  t2040 = t2035 + t2039;
  t2057 = -0.0265*t2034;
  t2061 = 0.2375*t2038;
  t2062 = t2057 + t2061;
  t2108 = -1.*t1836*t1855;
  t2115 = t1953*t1922;
  t2119 = t2108 + t2115;
  t2089 = t2018*t2051;
  t2042 = t2002*t2040;
  t2074 = t2051*t2062;
  t2077 = t2018*t2002;
  t2082 = t2051*t2038;
  t2084 = t2077 + t2082;
  t2085 = 0.0325*t2084;
  t2090 = -1.*t2002*t2038;
  t2093 = t2089 + t2090;
  t2094 = -0.0265*t2093;
  t2137 = -0.0265*t1855;
  t2142 = -0.0695*t1922;
  t2145 = t2137 + t2142;
  t2151 = 0.0695*t1855;
  t2157 = t2151 + t1925;
  t2107 = t2051*t2040;
  t2122 = t2119*t2062;
  t2123 = t2018*t2119;
  t2128 = -1.*t2051*t2038;
  t2129 = t2123 + t2128;
  t2132 = -0.0265*t2129;
  t2133 = t2119*t2038;
  t2134 = t2089 + t2133;
  t2135 = 0.0325*t2134;
  t2165 = t1855*t1953;
  t2166 = t1836*t1922;
  t2167 = t2165 + t2166;
  t2169 = -0.0265*t2018;
  t2172 = -0.2375*t2038;
  t2174 = t2169 + t2172;
  t2179 = 0.2375*t2018;
  t2184 = t2179 + t2039;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=1.;
  p_output1[4]=t1836*t1932 - 1.*t1953*t1977 + t2042 + t2074 + t2085 + t2094;
  p_output1[5]=-1.*t1932*t1953 - 1.*t1836*t1977 + t2107 + t2122 + t2132 + t2135;
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=t2042 + t2074 + t2085 + t2094 + t1953*t2145 + t1836*t2157;
  p_output1[11]=t2107 + t2122 + t2132 + t2135 + t1836*t2145 - 1.*t1953*t2157;
  p_output1[12]=-0.0265*(t2090 - 1.*t2018*t2167) + 0.0325*(t2077 - 1.*t2038*t2167) + t2167*t2174 + t2002*t2184;
  p_output1[13]=0.0325*t2093 - 0.0265*(-1.*t2002*t2018 + t2128) + t2002*t2174 + t2051*t2184;
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
