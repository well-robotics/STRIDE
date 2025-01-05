/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:32 GMT-05:00
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
  double t1864;
  double t1885;
  double t1890;
  double t1899;
  double t1858;
  double t1954;
  double t2005;
  double t2007;
  double t2009;
  double t2023;
  double t1994;
  double t1999;
  double t2002;
  double t2036;
  double t2037;
  double t2038;
  double t1893;
  double t1912;
  double t1922;
  double t1958;
  double t1973;
  double t1979;
  double t2013;
  double t2027;
  double t2033;
  double t2041;
  double t2058;
  double t2064;
  double t2088;
  double t2113;
  double t2114;
  double t2117;
  t1864 = Cos(var1[3]);
  t1885 = -1.*t1864;
  t1890 = 1. + t1885;
  t1899 = Sin(var1[3]);
  t1858 = Cos(var1[2]);
  t1954 = Sin(var1[2]);
  t2005 = Cos(var1[4]);
  t2007 = -1.*t2005;
  t2009 = 1. + t2007;
  t2023 = Sin(var1[4]);
  t1994 = t1858*t1864;
  t1999 = t1954*t1899;
  t2002 = t1994 + t1999;
  t2036 = t1864*t1954;
  t2037 = -1.*t1858*t1899;
  t2038 = t2036 + t2037;
  t1893 = -0.0265*t1890;
  t1912 = -0.0695*t1899;
  t1922 = t1893 + t1912;
  t1958 = -0.0695*t1890;
  t1973 = 0.0265*t1899;
  t1979 = t1958 + t1973;
  t2013 = -0.0265*t2009;
  t2027 = -0.2375*t2023;
  t2033 = t2013 + t2027;
  t2041 = -0.2375*t2009;
  t2058 = 0.0265*t2023;
  t2064 = t2041 + t2058;
  t2088 = t2005*t2002;
  t2113 = -1.*t1864*t1954;
  t2114 = t1858*t1899;
  t2117 = t2113 + t2114;
  p_output1[0]=t1858*t1922 + t1954*t1979 + t2002*t2033 - 0.0115*(-1.*t2002*t2023 + t2005*t2038) + t2038*t2064 - 0.0265*(t2023*t2038 + t2088) + var1[0];
  p_output1[1]=-1.*t1922*t1954 + t1858*t1979 + t2002*t2064 + t2033*t2117 - 0.0265*(t2002*t2023 + t2005*t2117) - 0.0115*(t2088 - 1.*t2023*t2117) + var1[1];
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
