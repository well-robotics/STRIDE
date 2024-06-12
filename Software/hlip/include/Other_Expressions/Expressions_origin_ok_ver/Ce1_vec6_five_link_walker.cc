/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:31 GMT-05:00
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
  double t859;
  double t892;
  double t830;
  double t850;
  double t829;
  double t849;
  double t860;
  double t866;
  double t882;
  double t883;
  double t884;
  double t916;
  double t917;
  double t918;
  double t836;
  double t855;
  double t856;
  double t858;
  double t894;
  double t895;
  double t911;
  double t926;
  double t927;
  double t928;
  double t919;
  double t922;
  double t923;
  double t924;
  double t925;
  double t929;
  double t935;
  double t936;
  double t937;
  double t912;
  double t913;
  double t914;
  double t938;
  double t939;
  double t940;
  double t941;
  double t942;
  double t943;
  double t930;
  double t961;
  double t962;
  double t963;
  double t964;
  double t931;
  double t965;
  double t950;
  double t951;
  double t952;
  double t915;
  double t920;
  double t974;
  double t946;
  double t947;
  double t948;
  double t975;
  double t976;
  double t977;
  double t994;
  double t995;
  double t996;
  double t986;
  double t987;
  double t988;
  double t990;
  double t991;
  double t992;
  double t993;
  double t997;
  double t1017;
  double t1018;
  double t1002;
  double t1020;
  double t1021;
  double t1004;
  t859 = Cos(var1[6]);
  t892 = Sin(var1[6]);
  t830 = Sin(var1[2]);
  t850 = Sin(var1[5]);
  t829 = Cos(var1[5]);
  t849 = Cos(var1[2]);
  t860 = -1.*t859;
  t866 = 1. + t860;
  t882 = 0.5*t866;
  t883 = 0.671885*t859;
  t884 = t882 + t883;
  t916 = t829*t859;
  t917 = -1.*t850*t892;
  t918 = t916 + t917;
  t836 = -1.*t829*t830;
  t855 = -1.*t849*t850;
  t856 = t836 + t855;
  t858 = 0.51185934*t856;
  t894 = t884*t892;
  t895 = -0.171885*t859*t892;
  t911 = t894 + t895;
  t926 = -1.*t859*t850;
  t927 = -1.*t829*t892;
  t928 = t926 + t927;
  t919 = t849*t918;
  t922 = t884*t859;
  t923 = Power(t892,2);
  t924 = 0.171885*t923;
  t925 = t922 + t924;
  t929 = t849*t928;
  t935 = t830*t928;
  t936 = t935 + t919;
  t937 = 0.85216*t911*t936;
  t912 = t859*t850;
  t913 = t829*t892;
  t914 = t912 + t913;
  t938 = -1.*t829*t859;
  t939 = t850*t892;
  t940 = t938 + t939;
  t941 = t830*t940;
  t942 = t929 + t941;
  t943 = 0.85216*t925*t942;
  t930 = -1.*t830*t918;
  t961 = -1.*t849*t829;
  t962 = t830*t850;
  t963 = t961 + t962;
  t964 = 0.51185934*t963;
  t931 = t929 + t930;
  t965 = -1.*t830*t928;
  t950 = Power(t859,2);
  t951 = -0.171885*t950;
  t952 = t922 + t951;
  t915 = -1.*t830*t914;
  t920 = t915 + t919;
  t974 = 0.85216*t911*t931;
  t946 = -1.*t884*t892;
  t947 = 0.171885*t859*t892;
  t948 = t946 + t947;
  t975 = t849*t940;
  t976 = t965 + t975;
  t977 = 0.85216*t925*t976;
  t994 = t829*t884;
  t995 = -0.171885*t850*t892;
  t996 = t994 + t995;
  t986 = -1.*t884*t850;
  t987 = -0.171885*t829*t892;
  t988 = t986 + t987;
  t990 = t884*t850;
  t991 = 0.171885*t829*t892;
  t992 = t990 + t991;
  t993 = t992*t918;
  t997 = t928*t996;
  t1017 = -0.171885*t859*t850;
  t1018 = t1017 + t987;
  t1002 = -1.*t928*t992;
  t1020 = 0.171885*t829*t859;
  t1021 = t1020 + t995;
  t1004 = -1.*t996*t940;
  p_output1[0]=var2[5]*(-0.5*(t858 + 0.85216*t911*t920 + 0.85216*t925*t931)*var2[2] - 0.5*(t858 + t937 + t943)*var2[5] - 0.5*(t937 + t943 + 0.85216*t936*t948 + 0.85216*(t849*t914 + t830*t918)*t952)*var2[6]);
  p_output1[1]=var2[5]*(-0.5*(0.85216*t911*(-1.*t849*t914 + t930) + t964 + 0.85216*t925*(-1.*t849*t918 + t965))*var2[2] - 0.5*(t964 + t974 + t977)*var2[5] - 0.5*(0.85216*t931*t948 + 0.85216*t920*t952 + t974 + t977)*var2[6]);
  p_output1[2]=var2[5]*(-0.5*(0.85216*t911*(t1002 + t1004 - 1.*t928*t988 - 1.*t918*t996) + 0.85216*t925*(t918*t988 + t993 + t914*t996 + t997))*var2[5] - 0.5*(0.85216*t911*(t1002 + t1004 - 1.*t1021*t918 - 1.*t1018*t928) + 0.85216*t948*(t914*t992 + t918*t996) + 0.85216*t952*(-1.*t918*t992 - 1.*t928*t996) + 0.85216*t925*(t1021*t914 + t1018*t918 + t993 + t997))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(1.70432*t925*t948 + 1.70432*t911*t952)*var2[5]*var2[6];
  p_output1[6]=-0.0732367608*t948*var2[5]*var2[6];
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

#include "Ce1_vec6_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
