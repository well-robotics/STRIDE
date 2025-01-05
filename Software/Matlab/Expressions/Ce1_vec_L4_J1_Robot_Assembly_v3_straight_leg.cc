/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:09 GMT-05:00
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
  double t330;
  double t290;
  double t323;
  double t340;
  double t366;
  double t368;
  double t369;
  double t326;
  double t348;
  double t354;
  double t370;
  double t371;
  double t372;
  double t373;
  double t374;
  double t376;
  double t403;
  double t404;
  double t409;
  double t410;
  double t411;
  double t412;
  double t415;
  double t416;
  double t417;
  double t418;
  double t425;
  double t427;
  double t436;
  double t443;
  double t444;
  double t445;
  double t429;
  double t430;
  double t431;
  double t432;
  double t435;
  double t448;
  double t453;
  double t454;
  double t459;
  double t466;
  double t468;
  double t472;
  double t488;
  double t489;
  double t481;
  double t482;
  double t483;
  double t505;
  double t507;
  double t508;
  t330 = Cos(var1[2]);
  t290 = Cos(var1[5]);
  t323 = Sin(var1[2]);
  t340 = Sin(var1[5]);
  t366 = t330*t290;
  t368 = -1.*t323*t340;
  t369 = t366 + t368;
  t326 = -1.*t290*t323;
  t348 = -1.*t330*t340;
  t354 = t326 + t348;
  t370 = 1.38102*t354*t369;
  t371 = t290*t323;
  t372 = t330*t340;
  t373 = t371 + t372;
  t374 = 1.38102*t373*t369;
  t376 = t370 + t374;
  t403 = Power(t354,2);
  t404 = 0.69051*t403;
  t409 = 0.69051*t354*t373;
  t410 = Power(t369,2);
  t411 = 0.69051*t410;
  t412 = -1.*t330*t290;
  t415 = t323*t340;
  t416 = t412 + t415;
  t417 = 0.69051*t369*t416;
  t418 = t404 + t409 + t411 + t417;
  t425 = -1.*t290;
  t427 = 1. + t425;
  t436 = -0.0265*t427;
  t443 = -0.025413*t290;
  t444 = -0.08282*t340;
  t445 = t436 + t443 + t444;
  t429 = -0.0695*t427;
  t430 = -0.15232*t290;
  t431 = -0.0010869999999999977*t340;
  t432 = t429 + t430 + t431;
  t435 = t290*t432;
  t448 = t445*t340;
  t453 = t435 + t448;
  t454 = 0.69051*t354*t453;
  t459 = -1.*t290*t445;
  t466 = t432*t340;
  t468 = t459 + t466;
  t472 = 0.69051*t369*t468;
  t488 = -0.08282*t290;
  t489 = t488 + t431;
  t481 = -0.0010869999999999977*t290;
  t482 = 0.08282*t340;
  t483 = t481 + t482;
  t505 = -0.0571880382*t354;
  t507 = -0.0007505843699999984*t369;
  t508 = t505 + t507;
  p_output1[0]=var2[0]*(-0.5*t376*var2[2] - 0.5*t376*var2[5]);
  p_output1[1]=var2[0]*(-0.5*t418*var2[2] - 0.5*t418*var2[5]);
  p_output1[2]=var2[0]*(-0.5*(t454 + t472)*var2[2] - 0.5*(t454 + t472 + 0.69051*t373*(t435 + t448 + t340*t483 - 1.*t290*t489) + 0.69051*t369*(-1.*t340*t432 + t290*t445 + t290*t483 + t340*t489))*var2[5]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(-0.5*t508*var2[2] - 0.5*t508*var2[5]);
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

#include "Ce1_vec_L4_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L4_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
