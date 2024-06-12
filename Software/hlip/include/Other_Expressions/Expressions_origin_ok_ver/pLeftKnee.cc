/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:46 GMT-05:00
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
  double t486;
  double t515;
  double t522;
  double t568;
  double t478;
  double t623;
  double t660;
  double t662;
  double t666;
  double t668;
  double t677;
  double t681;
  double t685;
  double t651;
  double t655;
  double t656;
  double t548;
  double t570;
  double t605;
  double t625;
  double t627;
  double t631;
  double t667;
  double t669;
  double t671;
  double t690;
  double t691;
  double t692;
  double t694;
  double t742;
  double t749;
  double t752;
  t486 = Cos(var1[5]);
  t515 = -1.*t486;
  t522 = 1. + t515;
  t568 = Sin(var1[5]);
  t478 = Sin(var1[2]);
  t623 = Cos(var1[2]);
  t660 = Cos(var1[6]);
  t662 = -1.*t660;
  t666 = 1. + t662;
  t668 = Sin(var1[6]);
  t677 = t623*t486;
  t681 = -1.*t478*t568;
  t685 = t677 + t681;
  t651 = t486*t478;
  t655 = t623*t568;
  t656 = t651 + t655;
  t548 = -0.0695*t522;
  t570 = -0.0265*t568;
  t605 = t548 + t570;
  t625 = -0.0265*t522;
  t627 = 0.0695*t568;
  t631 = t625 + t627;
  t667 = -0.2375*t666;
  t669 = -0.0265*t668;
  t671 = t667 + t669;
  t690 = -0.0265*t666;
  t691 = 0.2375*t668;
  t692 = t690 + t691;
  t694 = t660*t685;
  t742 = -1.*t486*t478;
  t749 = -1.*t623*t568;
  t752 = t742 + t749;
  p_output1[0]=t478*t605 + t623*t631 + t656*t671 - 0.2375*(t656*t660 + t668*t685) + t685*t692 - 0.0265*(-1.*t656*t668 + t694) + var1[0];
  p_output1[1]=0.0293;
  p_output1[2]=t605*t623 - 1.*t478*t631 + t671*t685 + t692*t752 - 0.0265*(-1.*t668*t685 + t660*t752) - 0.2375*(t694 + t668*t752) + var1[1];
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

#include "pLeftKnee.hh"

namespace SymFunction
{

void pLeftKnee_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
