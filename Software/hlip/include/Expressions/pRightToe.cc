/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:41 GMT-05:00
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
  double t1807;
  double t1828;
  double t1833;
  double t1842;
  double t1801;
  double t1897;
  double t1947;
  double t1950;
  double t1952;
  double t1966;
  double t1932;
  double t1938;
  double t1943;
  double t1979;
  double t1980;
  double t1981;
  double t1836;
  double t1855;
  double t1865;
  double t1901;
  double t1916;
  double t1922;
  double t1953;
  double t1970;
  double t1976;
  double t1982;
  double t2000;
  double t2002;
  double t2031;
  double t2056;
  double t2057;
  double t2060;
  t1807 = Cos(var1[3]);
  t1828 = -1.*t1807;
  t1833 = 1. + t1828;
  t1842 = Sin(var1[3]);
  t1801 = Cos(var1[2]);
  t1897 = Sin(var1[2]);
  t1947 = Cos(var1[4]);
  t1950 = -1.*t1947;
  t1952 = 1. + t1950;
  t1966 = Sin(var1[4]);
  t1932 = t1801*t1807;
  t1938 = t1897*t1842;
  t1943 = t1932 + t1938;
  t1979 = t1807*t1897;
  t1980 = -1.*t1801*t1842;
  t1981 = t1979 + t1980;
  t1836 = -0.0265*t1833;
  t1855 = -0.0695*t1842;
  t1865 = t1836 + t1855;
  t1901 = -0.0695*t1833;
  t1916 = 0.0265*t1842;
  t1922 = t1901 + t1916;
  t1953 = -0.0265*t1952;
  t1970 = -0.2375*t1966;
  t1976 = t1953 + t1970;
  t1982 = -0.2375*t1952;
  t2000 = 0.0265*t1966;
  t2002 = t1982 + t2000;
  t2031 = t1947*t1943;
  t2056 = -1.*t1807*t1897;
  t2057 = t1801*t1842;
  t2060 = t2056 + t2057;
  p_output1[0]=t1801*t1865 + t1897*t1922 + t1943*t1976 + 0.0325*(-1.*t1943*t1966 + t1947*t1981) + t1981*t2002 - 0.0265*(t1966*t1981 + t2031) + var1[0];
  p_output1[1]=-1.*t1865*t1897 + t1801*t1922 + t1943*t2002 + t1976*t2060 - 0.0265*(t1943*t1966 + t1947*t2060) + 0.0325*(t2031 - 1.*t1966*t2060) + var1[1];
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

#include "pRightToe.hh"

namespace SymFunction
{

void pRightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
