/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:40 GMT-05:00
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
  double t1775;
  double t1789;
  double t1797;
  double t1803;
  double t1761;
  double t1834;
  double t1885;
  double t1888;
  double t1894;
  double t1898;
  double t1925;
  double t1926;
  double t1927;
  double t1865;
  double t1869;
  double t1884;
  double t1801;
  double t1807;
  double t1828;
  double t1836;
  double t1842;
  double t1851;
  double t1897;
  double t1901;
  double t1916;
  double t1932;
  double t1938;
  double t1943;
  double t1949;
  double t1985;
  double t1988;
  double t1995;
  t1775 = Cos(var1[5]);
  t1789 = -1.*t1775;
  t1797 = 1. + t1789;
  t1803 = Sin(var1[5]);
  t1761 = Sin(var1[2]);
  t1834 = Cos(var1[2]);
  t1885 = Cos(var1[6]);
  t1888 = -1.*t1885;
  t1894 = 1. + t1888;
  t1898 = Sin(var1[6]);
  t1925 = t1834*t1775;
  t1926 = -1.*t1761*t1803;
  t1927 = t1925 + t1926;
  t1865 = t1775*t1761;
  t1869 = t1834*t1803;
  t1884 = t1865 + t1869;
  t1801 = -0.0695*t1797;
  t1807 = -0.0265*t1803;
  t1828 = t1801 + t1807;
  t1836 = -0.0265*t1797;
  t1842 = 0.0695*t1803;
  t1851 = t1836 + t1842;
  t1897 = -0.2375*t1894;
  t1901 = -0.0265*t1898;
  t1916 = t1897 + t1901;
  t1932 = -0.0265*t1894;
  t1938 = 0.2375*t1898;
  t1943 = t1932 + t1938;
  t1949 = t1885*t1927;
  t1985 = -1.*t1775*t1761;
  t1988 = -1.*t1834*t1803;
  t1995 = t1985 + t1988;
  p_output1[0]=t1761*t1828 + t1834*t1851 + t1884*t1916 + 0.0325*(t1884*t1885 + t1898*t1927) + t1927*t1943 - 0.0265*(-1.*t1884*t1898 + t1949) + var1[0];
  p_output1[1]=t1828*t1834 - 1.*t1761*t1851 + t1916*t1927 + t1943*t1995 - 0.0265*(-1.*t1898*t1927 + t1885*t1995) + 0.0325*(t1949 + t1898*t1995) + var1[1];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "pLeftToe.hh"

namespace SymFunction
{

void pLeftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
