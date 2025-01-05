/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:12 GMT-05:00
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
  double t354;
  double t366;
  double t369;
  double t370;
  double t371;
  double t372;
  double t373;
  double t404;
  double t409;
  double t410;
  double t374;
  double t383;
  double t391;
  double t394;
  double t398;
  double t403;
  double t411;
  double t415;
  double t416;
  double t417;
  double t419;
  double t421;
  double t422;
  double t436;
  double t443;
  double t444;
  double t479;
  double t480;
  double t490;
  double t491;
  double t492;
  double t493;
  double t481;
  double t482;
  double t486;
  double t487;
  double t488;
  double t494;
  double t495;
  double t496;
  double t497;
  double t498;
  double t500;
  double t501;
  double t517;
  double t518;
  double t512;
  double t513;
  double t514;
  double t530;
  double t531;
  double t532;
  t354 = Cos(var1[5]);
  t366 = Sin(var1[2]);
  t369 = -1.*t354*t366;
  t370 = Cos(var1[2]);
  t371 = Sin(var1[5]);
  t372 = -1.*t370*t371;
  t373 = t369 + t372;
  t404 = t370*t354;
  t409 = -1.*t366*t371;
  t410 = t404 + t409;
  t374 = Power(t373,2);
  t383 = 0.69051*t374;
  t391 = t354*t366;
  t394 = t370*t371;
  t398 = t391 + t394;
  t403 = 0.69051*t373*t398;
  t411 = Power(t410,2);
  t415 = 0.69051*t411;
  t416 = -1.*t370*t354;
  t417 = t366*t371;
  t419 = t416 + t417;
  t421 = 0.69051*t410*t419;
  t422 = t383 + t403 + t415 + t421;
  t436 = 1.38102*t373*t410;
  t443 = 1.38102*t373*t419;
  t444 = t436 + t443;
  t479 = -1.*t354;
  t480 = 1. + t479;
  t490 = -0.0265*t480;
  t491 = -0.025413*t354;
  t492 = -0.08282*t371;
  t493 = t490 + t491 + t492;
  t481 = -0.0695*t480;
  t482 = -0.15232*t354;
  t486 = -0.0010869999999999977*t371;
  t487 = t481 + t482 + t486;
  t488 = t354*t487;
  t494 = t493*t371;
  t495 = t488 + t494;
  t496 = 0.69051*t419*t495;
  t497 = -1.*t354*t493;
  t498 = t487*t371;
  t500 = t497 + t498;
  t501 = 0.69051*t373*t500;
  t517 = -0.08282*t354;
  t518 = t517 + t486;
  t512 = -0.0010869999999999977*t354;
  t513 = 0.08282*t371;
  t514 = t512 + t513;
  t530 = -0.0007505843699999984*t373;
  t531 = -0.0571880382*t419;
  t532 = t530 + t531;
  p_output1[0]=var2[1]*(-0.5*t422*var2[2] - 0.5*t422*var2[5]);
  p_output1[1]=var2[1]*(-0.5*t444*var2[2] - 0.5*t444*var2[5]);
  p_output1[2]=var2[1]*(-0.5*(t496 + t501)*var2[2] - 0.5*(t496 + t501 + 0.69051*t410*(t488 + t494 + t371*t514 - 1.*t354*t518) + 0.69051*t373*(-1.*t371*t487 + t354*t493 + t354*t514 + t371*t518))*var2[5]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(-0.5*t532*var2[2] - 0.5*t532*var2[5]);
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

#include "Ce1_vec_L4_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L4_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
