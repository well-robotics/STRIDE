/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:48 GMT-05:00
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
  double t570;
  double t605;
  double t606;
  double t627;
  double t548;
  double t667;
  double t693;
  double t695;
  double t699;
  double t710;
  double t690;
  double t691;
  double t692;
  double t732;
  double t736;
  double t739;
  double t625;
  double t641;
  double t651;
  double t669;
  double t671;
  double t672;
  double t701;
  double t712;
  double t721;
  double t740;
  double t759;
  double t760;
  double t768;
  double t813;
  double t814;
  double t815;
  t570 = Cos(var1[3]);
  t605 = -1.*t570;
  t606 = 1. + t605;
  t627 = Sin(var1[3]);
  t548 = Cos(var1[2]);
  t667 = Sin(var1[2]);
  t693 = Cos(var1[4]);
  t695 = -1.*t693;
  t699 = 1. + t695;
  t710 = Sin(var1[4]);
  t690 = t548*t570;
  t691 = t667*t627;
  t692 = t690 + t691;
  t732 = t570*t667;
  t736 = -1.*t548*t627;
  t739 = t732 + t736;
  t625 = -0.0265*t606;
  t641 = -0.0695*t627;
  t651 = t625 + t641;
  t669 = -0.0695*t606;
  t671 = 0.0265*t627;
  t672 = t669 + t671;
  t701 = -0.0265*t699;
  t712 = -0.2375*t710;
  t721 = t701 + t712;
  t740 = -0.2375*t699;
  t759 = 0.0265*t710;
  t760 = t740 + t759;
  t768 = t693*t692;
  t813 = -1.*t570*t667;
  t814 = t548*t627;
  t815 = t813 + t814;
  p_output1[0]=t548*t651 + t667*t672 + t692*t721 - 0.2375*(-1.*t692*t710 + t693*t739) + t739*t760 - 0.0265*(t710*t739 + t768) + var1[0];
  p_output1[1]=-0.0293;
  p_output1[2]=-1.*t651*t667 + t548*t672 + t692*t760 + t721*t815 - 0.0265*(t692*t710 + t693*t815) - 0.2375*(t768 - 1.*t710*t815) + var1[1];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 3, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "pRightKnee.hh"

namespace SymFunction
{

void pRightKnee_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
