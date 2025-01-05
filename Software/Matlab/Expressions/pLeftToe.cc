/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:30 GMT-05:00
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
  double t1832;
  double t1846;
  double t1854;
  double t1861;
  double t1825;
  double t1891;
  double t1942;
  double t1950;
  double t1951;
  double t1955;
  double t1982;
  double t1983;
  double t1984;
  double t1922;
  double t1939;
  double t1941;
  double t1858;
  double t1864;
  double t1885;
  double t1893;
  double t1899;
  double t1908;
  double t1954;
  double t1958;
  double t1973;
  double t1994;
  double t1999;
  double t2002;
  double t2006;
  double t2042;
  double t2045;
  double t2052;
  t1832 = Cos(var1[5]);
  t1846 = -1.*t1832;
  t1854 = 1. + t1846;
  t1861 = Sin(var1[5]);
  t1825 = Sin(var1[2]);
  t1891 = Cos(var1[2]);
  t1942 = Cos(var1[6]);
  t1950 = -1.*t1942;
  t1951 = 1. + t1950;
  t1955 = Sin(var1[6]);
  t1982 = t1891*t1832;
  t1983 = -1.*t1825*t1861;
  t1984 = t1982 + t1983;
  t1922 = t1832*t1825;
  t1939 = t1891*t1861;
  t1941 = t1922 + t1939;
  t1858 = -0.0695*t1854;
  t1864 = -0.0265*t1861;
  t1885 = t1858 + t1864;
  t1893 = -0.0265*t1854;
  t1899 = 0.0695*t1861;
  t1908 = t1893 + t1899;
  t1954 = -0.2375*t1951;
  t1958 = -0.0265*t1955;
  t1973 = t1954 + t1958;
  t1994 = -0.0265*t1951;
  t1999 = 0.2375*t1955;
  t2002 = t1994 + t1999;
  t2006 = t1942*t1984;
  t2042 = -1.*t1832*t1825;
  t2045 = -1.*t1891*t1861;
  t2052 = t2042 + t2045;
  p_output1[0]=t1825*t1885 + t1891*t1908 + t1941*t1973 - 0.0115*(t1941*t1942 + t1955*t1984) + t1984*t2002 - 0.0265*(-1.*t1941*t1955 + t2006) + var1[0];
  p_output1[1]=t1885*t1891 - 1.*t1825*t1908 + t1973*t1984 + t2002*t2052 - 0.0265*(-1.*t1955*t1984 + t1942*t2052) - 0.0115*(t2006 + t1955*t2052) + var1[1];
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
