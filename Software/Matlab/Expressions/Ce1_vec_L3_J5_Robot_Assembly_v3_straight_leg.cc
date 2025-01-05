/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:02 GMT-05:00
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
  double t269;
  double t257;
  double t260;
  double t270;
  double t287;
  double t251;
  double t290;
  double t318;
  double t321;
  double t262;
  double t271;
  double t278;
  double t371;
  double t372;
  double t373;
  double t374;
  double t376;
  double t379;
  double t383;
  double t387;
  double t391;
  double t394;
  double t398;
  double t402;
  double t340;
  double t348;
  double t349;
  double t359;
  double t282;
  double t323;
  double t326;
  double t423;
  double t424;
  double t425;
  double t427;
  double t428;
  double t437;
  double t441;
  double t436;
  double t443;
  double t444;
  double t445;
  double t446;
  double t447;
  double t448;
  double t453;
  double t454;
  double t458;
  double t459;
  double t433;
  double t434;
  double t472;
  double t473;
  double t474;
  double t475;
  double t476;
  double t480;
  double t481;
  double t482;
  double t483;
  double t484;
  double t466;
  double t467;
  double t435;
  double t461;
  double t464;
  double t485;
  double t496;
  double t497;
  double t498;
  double t501;
  double t504;
  double t489;
  double t508;
  double t509;
  double t510;
  double t500;
  double t505;
  double t506;
  double t491;
  t269 = Cos(var1[3]);
  t257 = Cos(var1[4]);
  t260 = Sin(var1[3]);
  t270 = Sin(var1[4]);
  t287 = Cos(var1[2]);
  t251 = Sin(var1[2]);
  t290 = t269*t257;
  t318 = -1.*t260*t270;
  t321 = t290 + t318;
  t262 = -1.*t257*t260;
  t271 = -1.*t269*t270;
  t278 = t262 + t271;
  t371 = t287*t278;
  t372 = t251*t321;
  t373 = t371 + t372;
  t374 = 0.015375074960000006*t373;
  t376 = t251*t278;
  t379 = -1.*t269*t257;
  t383 = t260*t270;
  t387 = t379 + t383;
  t391 = t287*t387;
  t394 = t376 + t391;
  t398 = 0.0002537424399999996*t394;
  t402 = t374 + t398;
  t340 = t257*t260;
  t348 = t269*t270;
  t349 = t340 + t348;
  t359 = -1.*t251*t321;
  t282 = -1.*t251*t278;
  t323 = t287*t321;
  t326 = t282 + t323;
  t423 = 0.015375074960000006*t326;
  t424 = -1.*t251*t387;
  t425 = t371 + t424;
  t427 = 0.0002537424399999996*t425;
  t428 = t423 + t427;
  t437 = -1.*t257;
  t441 = 1. + t437;
  t436 = -0.0695*t260;
  t443 = -0.2375*t441;
  t444 = -0.314514*t257;
  t445 = 0.0012709999999999978*t270;
  t446 = t443 + t444 + t445;
  t447 = -1.*t260*t446;
  t448 = -0.0265*t441;
  t453 = -0.025229*t257;
  t454 = 0.07701400000000003*t270;
  t458 = t448 + t453 + t454;
  t459 = t269*t458;
  t433 = -1.*t269;
  t434 = 1. + t433;
  t472 = -0.0695*t269;
  t473 = -0.0265*t260;
  t474 = -1.*t269*t446;
  t475 = -1.*t260*t458;
  t476 = t472 + t473 + t474 + t475;
  t480 = -0.0695*t434;
  t481 = 0.0265*t260;
  t482 = t269*t446;
  t483 = t260*t458;
  t484 = t480 + t481 + t482 + t483;
  t466 = 0.0265*t269;
  t467 = t466 + t436 + t447 + t459;
  t435 = -0.0265*t434;
  t461 = t435 + t436 + t447 + t459;
  t464 = -1.*t461*t278;
  t485 = -1.*t484*t321;
  t496 = 0.07701400000000003*t257;
  t497 = -0.0012709999999999978*t270;
  t498 = t496 + t497;
  t501 = 0.0012709999999999978*t257;
  t504 = t501 + t454;
  t489 = t484*t278;
  t508 = t269*t498;
  t509 = -1.*t260*t504;
  t510 = t508 + t509;
  t500 = t260*t498;
  t505 = t269*t504;
  t506 = t500 + t505;
  t491 = t461*t387;
  p_output1[0]=var2[4]*(-0.5*(0.0002537424399999996*t326 + 0.015375074960000006*(t287*t349 + t359))*var2[2] - 0.5*t402*var2[3] - 0.5*t402*var2[4]);
  p_output1[1]=var2[4]*(-0.5*(0.015375074960000006*(-1.*t287*t321 - 1.*t251*t349) + 0.0002537424399999996*(-1.*t278*t287 + t359))*var2[2] - 0.5*t428*var2[3] - 0.5*t428*var2[4]);
  p_output1[2]=var2[4]*(-0.5*(0.0002537424399999996*(t464 - 1.*t349*t467 - 1.*t321*t476 + t485) + 0.015375074960000006*(t321*t467 + t278*t476 + t489 + t491))*var2[3] - 0.5*(0.015375074960000006*(t489 + t491 + t321*t506 + t278*t510) + 0.0002537424399999996*(t464 + t485 - 1.*t349*t506 - 1.*t321*t510))*var2[4]);
  p_output1[3]=-0.5*(0.015375074960000006*(0.0265*t257 + 0.0695*t270 + t270*t446 + t257*t458 + t270*t498 - 1.*t257*t504) + 0.0002537424399999996*(0.0695*t257 - 0.0265*t270 + t257*t446 - 1.*t270*t458 + t257*t498 + t270*t504))*Power(var2[4],2);
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

#include "Ce1_vec_L3_J5_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L3_J5_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
