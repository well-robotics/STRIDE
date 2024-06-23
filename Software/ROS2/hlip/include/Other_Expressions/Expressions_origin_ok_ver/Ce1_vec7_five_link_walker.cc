/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:32 GMT-05:00
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
  double t944;
  double t932;
  double t933;
  double t945;
  double t921;
  double t934;
  double t949;
  double t953;
  double t954;
  double t955;
  double t967;
  double t968;
  double t969;
  double t970;
  double t971;
  double t956;
  double t957;
  double t958;
  double t980;
  double t984;
  double t985;
  double t1001;
  double t1003;
  double t1005;
  double t1006;
  double t1007;
  double t1016;
  double t1019;
  double t1022;
  double t1009;
  double t1012;
  double t1013;
  double t1014;
  double t1015;
  double t1023;
  double t1024;
  double t1025;
  double t1026;
  t944 = Cos(var1[5]);
  t932 = Cos(var1[6]);
  t933 = Sin(var1[5]);
  t945 = Sin(var1[6]);
  t921 = Cos(var1[2]);
  t934 = -1.*t932*t933;
  t949 = -1.*t944*t945;
  t953 = t934 + t949;
  t954 = t921*t953;
  t955 = Sin(var1[2]);
  t967 = -1.*t944*t932;
  t968 = t933*t945;
  t969 = t967 + t968;
  t970 = t955*t969;
  t971 = t954 + t970;
  t956 = t944*t932;
  t957 = -1.*t933*t945;
  t958 = t956 + t957;
  t980 = -1.*t955*t953;
  t984 = t921*t969;
  t985 = t980 + t984;
  t1001 = -1.*t932;
  t1003 = 1. + t1001;
  t1005 = 0.5*t1003;
  t1006 = 0.671885*t932;
  t1007 = t1005 + t1006;
  t1016 = t944*t1007;
  t1019 = -0.171885*t933*t945;
  t1022 = t1016 + t1019;
  t1009 = -0.171885*t944*t945;
  t1012 = t1007*t933;
  t1013 = 0.171885*t944*t945;
  t1014 = t1012 + t1013;
  t1015 = t1014*t958;
  t1023 = t953*t1022;
  t1024 = t932*t933;
  t1025 = t944*t945;
  t1026 = t1024 + t1025;
  p_output1[0]=var2[6]*(-0.0732367608*(t954 - 1.*t955*t958)*var2[2] - 0.0732367608*t971*var2[5] - 0.0732367608*t971*var2[6]);
  p_output1[1]=var2[6]*(-0.0732367608*(-1.*t921*t958 + t980)*var2[2] - 0.0732367608*t985*var2[5] - 0.0732367608*t985*var2[6]);
  p_output1[2]=var2[6]*(-0.0732367608*(t1015 + t1023 + t1022*t1026 + (t1009 - 1.*t1007*t933)*t958)*var2[5] - 0.0732367608*(t1015 + t1023 + t1026*(t1019 + 0.171885*t932*t944) + (t1009 - 0.171885*t932*t933)*t958)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.0732367608*(-1.*t1007*t945 + 0.171885*t932*t945)*Power(var2[6],2);
  p_output1[6]=0;
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

#include "Ce1_vec7_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
