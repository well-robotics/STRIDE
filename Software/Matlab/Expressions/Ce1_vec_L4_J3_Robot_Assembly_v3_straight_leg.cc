/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:14 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t373;
  double t394;
  double t404;
  double t410;
  double t391;
  double t374;
  double t423;
  double t429;
  double t430;
  double t431;
  double t411;
  double t415;
  double t417;
  double t419;
  double t383;
  double t398;
  double t403;
  double t421;
  double t436;
  double t443;
  double t453;
  double t466;
  double t468;
  double t473;
  double t481;
  double t482;
  double t486;
  double t490;
  double t512;
  double t513;
  double t498;
  double t500;
  double t505;
  double t530;
  double t531;
  double t533;
  double t534;
  double t535;
  double t495;
  double t507;
  double t511;
  double t515;
  double t516;
  double t522;
  double t523;
  double t524;
  t373 = Cos(var1[5]);
  t394 = Sin(var1[5]);
  t404 = -1.*t373;
  t410 = 1. + t404;
  t391 = Cos(var1[2]);
  t374 = Sin(var1[2]);
  t423 = -0.0265*t410;
  t429 = -0.025413*t373;
  t430 = -0.08282*t394;
  t431 = t423 + t429 + t430;
  t411 = -0.0695*t410;
  t415 = -0.15232*t373;
  t417 = -0.0010869999999999977*t394;
  t419 = t411 + t415 + t417;
  t383 = -1.*t373*t374;
  t398 = -1.*t391*t394;
  t403 = t383 + t398;
  t421 = t373*t419;
  t436 = t431*t394;
  t443 = t421 + t436;
  t453 = 0.69051*t403*t443;
  t466 = t391*t373;
  t468 = -1.*t374*t394;
  t473 = t466 + t468;
  t481 = -1.*t373*t431;
  t482 = t419*t394;
  t486 = t481 + t482;
  t490 = 0.69051*t473*t486;
  t512 = -0.08282*t373;
  t513 = t512 + t417;
  t498 = -0.0010869999999999977*t373;
  t500 = 0.08282*t394;
  t505 = t498 + t500;
  t530 = -1.*t391*t373;
  t531 = t374*t394;
  t533 = t530 + t531;
  t534 = 0.69051*t533*t443;
  t535 = 0.69051*t403*t486;
  t495 = t373*t431;
  t507 = t373*t505;
  t511 = -1.*t419*t394;
  t515 = t513*t394;
  t516 = t495 + t507 + t511 + t515;
  t522 = -1.*t373*t513;
  t523 = t505*t394;
  t524 = t421 + t522 + t436 + t523;
  p_output1[0]=var2[2]*(-0.5*(t453 + t490)*var2[2] - 0.5*(t453 + t490 + 0.69051*t473*t516 + 0.69051*(t373*t374 + t391*t394)*t524)*var2[5]);
  p_output1[1]=var2[2]*(-0.5*(t534 + t535)*var2[2] - 0.5*(0.69051*t403*t516 + 0.69051*t473*t524 + t534 + t535)*var2[5]);
  p_output1[2]=-0.5*(1.38102*t443*t516 + 1.38102*t486*t524)*var2[2]*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(-0.0571880382*t516 - 0.0007505843699999984*t524)*var2[2]*var2[5];
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

#include "Ce1_vec_L4_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L4_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
