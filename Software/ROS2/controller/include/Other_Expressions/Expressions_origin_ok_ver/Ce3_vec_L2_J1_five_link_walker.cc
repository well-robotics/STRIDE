/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:20 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t854;
  double t843;
  double t844;
  double t861;
  double t846;
  double t863;
  double t869;
  double t879;
  double t881;
  double t883;
  double t910;
  double t911;
  double t916;
  double t960;
  double t961;
  double t963;
  double t951;
  double t952;
  double t953;
  double t877;
  double t889;
  double t891;
  double t893;
  double t909;
  double t918;
  double t920;
  double t922;
  double t923;
  double t925;
  double t926;
  double t928;
  double t935;
  double t940;
  double t941;
  double t942;
  double t943;
  double t949;
  double t950;
  double t955;
  double t964;
  double t965;
  double t967;
  double t969;
  double t970;
  double t971;
  double t973;
  double t982;
  double t984;
  double t985;
  t854 = Cos(var1[2]);
  t843 = Cos(var1[3]);
  t844 = Sin(var1[2]);
  t861 = Sin(var1[3]);
  t846 = -1.*t843*t844;
  t863 = -1.*t854*t861;
  t869 = t846 + t863;
  t879 = t854*t843;
  t881 = -1.*t844*t861;
  t883 = t879 + t881;
  t910 = t843*t844;
  t911 = t854*t861;
  t916 = t910 + t911;
  t960 = -0.00102*t843;
  t961 = -0.078722*t861;
  t963 = t960 + t961;
  t951 = -0.078722*t843;
  t952 = 0.00102*t861;
  t953 = t951 + t952;
  t877 = -0.05089692188*t869;
  t889 = 0.0006594708000000001*t883;
  t891 = t877 + t889;
  t893 = 0.5*var2[3]*t891;
  t909 = 1.29308*t869*t883;
  t918 = 1.29308*t916*t883;
  t920 = t909 + t918;
  t922 = 0.5*var2[0]*t920;
  t923 = Power(t869,2);
  t925 = 0.64654*t923;
  t926 = 0.64654*t869*t916;
  t928 = Power(t883,2);
  t935 = 0.64654*t928;
  t940 = -1.*t854*t843;
  t941 = t844*t861;
  t942 = t940 + t941;
  t943 = 0.64654*t883*t942;
  t949 = t925 + t926 + t935 + t943;
  t950 = 0.5*var2[1]*t949;
  t955 = t843*t953;
  t964 = t963*t861;
  t965 = t955 + t964;
  t967 = 0.64654*t869*t965;
  t969 = -1.*t843*t963;
  t970 = t953*t861;
  t971 = t969 + t970;
  t973 = 0.64654*t883*t971;
  t982 = 0.00102*t843;
  t984 = 0.078722*t861;
  t985 = t982 + t984;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(t893 + t922 + t950 + 0.5*(t967 + t973)*var2[2]);
  p_output1[3]=var2[0]*(t893 + t922 + t950 + 0.5*(t967 + t973 + 0.64654*t883*(t843*t963 + t843*t985) + 0.64654*t916*(t964 + t861*t985))*var2[2]);
  p_output1[4]=0;
  p_output1[5]=0;
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

#include "Ce3_vec_L2_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L2_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
