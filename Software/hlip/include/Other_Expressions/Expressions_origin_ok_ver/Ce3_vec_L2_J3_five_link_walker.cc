/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:24 GMT-05:00
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
  double t911;
  double t925;
  double t923;
  double t910;
  double t961;
  double t964;
  double t965;
  double t943;
  double t951;
  double t955;
  double t976;
  double t981;
  double t982;
  double t960;
  double t970;
  double t971;
  double t984;
  double t985;
  double t987;
  double t1015;
  double t1020;
  double t1021;
  double t1014;
  double t1023;
  double t1024;
  double t1002;
  double t1005;
  double t1006;
  double t1001;
  double t1009;
  double t1026;
  double t1027;
  double t918;
  double t935;
  double t942;
  double t975;
  double t990;
  t911 = Cos(var1[3]);
  t925 = Sin(var1[3]);
  t923 = Sin(var1[2]);
  t910 = Cos(var1[2]);
  t961 = -0.00102*t911;
  t964 = -0.078722*t925;
  t965 = t961 + t964;
  t943 = -0.078722*t911;
  t951 = 0.00102*t925;
  t955 = t943 + t951;
  t976 = -1.*t911*t923;
  t981 = -1.*t910*t925;
  t982 = t976 + t981;
  t960 = t911*t955;
  t970 = t965*t925;
  t971 = t960 + t970;
  t984 = -1.*t911*t965;
  t985 = t955*t925;
  t987 = t984 + t985;
  t1015 = 0.00102*t911;
  t1020 = 0.078722*t925;
  t1021 = t1015 + t1020;
  t1014 = t911*t965;
  t1023 = t911*t1021;
  t1024 = t1014 + t1023;
  t1002 = t910*t911;
  t1005 = -1.*t923*t925;
  t1006 = t1002 + t1005;
  t1001 = 0.64654*t982*t971;
  t1009 = 0.64654*t1006*t987;
  t1026 = t1021*t925;
  t1027 = t970 + t1026;
  t918 = -1.*t910*t911;
  t935 = t923*t925;
  t942 = t918 + t935;
  t975 = 0.64654*t942*t971;
  t990 = 0.64654*t982*t987;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(t1001 + t1009)*var2[0] + 0.5*(t975 + t990)*var2[1])*var2[2];
  p_output1[3]=var2[2]*(0.5*(t1001 + t1009 + 0.64654*t1006*t1024 + 0.64654*t1027*(t911*t923 + t910*t925))*var2[0] + 0.5*(0.64654*t1006*t1027 + t975 + 0.64654*t1024*t982 + t990)*var2[1] + 0.5*(1.29308*t1024*t971 + 1.29308*t1027*t987)*var2[2] + 0.5*(-0.05089692188*t1024 + 0.0006594708000000001*t1027)*var2[3]);
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

#include "Ce3_vec_L2_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L2_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
