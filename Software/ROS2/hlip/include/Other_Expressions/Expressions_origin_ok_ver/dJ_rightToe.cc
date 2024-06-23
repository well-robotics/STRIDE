/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:34 GMT-05:00
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
  double t327;
  double t341;
  double t342;
  double t356;
  double t267;
  double t367;
  double t380;
  double t381;
  double t382;
  double t399;
  double t373;
  double t376;
  double t378;
  double t411;
  double t414;
  double t415;
  double t360;
  double t393;
  double t403;
  double t407;
  double t419;
  double t420;
  double t423;
  double t457;
  double t458;
  double t461;
  double t463;
  double t464;
  double t465;
  double t467;
  double t472;
  double t477;
  double t478;
  double t480;
  double t481;
  double t482;
  double t484;
  double t434;
  double t435;
  double t436;
  double t349;
  double t363;
  double t368;
  double t369;
  double t371;
  double t428;
  double t471;
  double t443;
  double t444;
  double t450;
  double t451;
  double t454;
  double t485;
  double t505;
  double t506;
  double t409;
  double t426;
  double t429;
  double t431;
  double t433;
  double t438;
  double t445;
  double t455;
  double t462;
  double t466;
  double t468;
  double t469;
  double t470;
  double t473;
  double t474;
  double t475;
  double t533;
  double t534;
  double t535;
  double t503;
  double t504;
  double t507;
  double t508;
  double t509;
  double t510;
  double t513;
  double t514;
  double t515;
  double t516;
  double t517;
  double t518;
  double t519;
  double t520;
  double t521;
  double t524;
  double t525;
  double t526;
  double t527;
  double t528;
  double t529;
  double t530;
  double t479;
  double t483;
  double t486;
  double t487;
  double t488;
  double t489;
  double t490;
  double t547;
  double t548;
  double t493;
  double t494;
  double t495;
  double t496;
  double t497;
  double t498;
  double t499;
  double t545;
  double t546;
  double t549;
  double t550;
  double t551;
  double t557;
  double t558;
  double t559;
  t327 = Cos(var1[3]);
  t341 = -1.*t327;
  t342 = 1. + t341;
  t356 = Sin(var1[3]);
  t267 = Cos(var1[2]);
  t367 = Sin(var1[2]);
  t380 = Cos(var1[4]);
  t381 = -1.*t380;
  t382 = 1. + t381;
  t399 = Sin(var1[4]);
  t373 = -1.*t267*t327;
  t376 = -1.*t367*t356;
  t378 = t373 + t376;
  t411 = -1.*t327*t367;
  t414 = t267*t356;
  t415 = t411 + t414;
  t360 = -0.0695*t356;
  t393 = -0.0265*t382;
  t403 = -0.2375*t399;
  t407 = t393 + t403;
  t419 = -0.2375*t382;
  t420 = 0.0265*t399;
  t423 = t419 + t420;
  t457 = t267*t327;
  t458 = t367*t356;
  t461 = t457 + t458;
  t463 = t327*t367;
  t464 = -1.*t267*t356;
  t465 = t463 + t464;
  t467 = t380*t461;
  t472 = -1.*t461*t399;
  t477 = 0.0265*t380;
  t478 = t477 + t403;
  t480 = -0.2375*t380;
  t481 = -0.0265*t399;
  t482 = t480 + t481;
  t484 = -1.*t415*t399;
  t434 = t380*t415;
  t435 = -1.*t378*t399;
  t436 = t434 + t435;
  t349 = -0.0265*t342;
  t363 = t349 + t360;
  t368 = -0.0695*t342;
  t369 = 0.0265*t356;
  t371 = t368 + t369;
  t428 = t380*t378;
  t471 = t380*t465;
  t443 = 0.0265*t327;
  t444 = t443 + t360;
  t450 = -0.0695*t327;
  t451 = -0.0265*t356;
  t454 = t450 + t451;
  t485 = t467 + t484;
  t505 = -1.*t465*t399;
  t506 = t428 + t505;
  t409 = t378*t407;
  t426 = t415*t423;
  t429 = t415*t399;
  t431 = t428 + t429;
  t433 = -0.0265*t431;
  t438 = 0.0225*t436;
  t445 = t267*t444;
  t455 = -1.*t367*t454;
  t462 = t461*t407;
  t466 = t465*t423;
  t468 = t465*t399;
  t469 = t467 + t468;
  t470 = -0.0265*t469;
  t473 = t471 + t472;
  t474 = 0.0225*t473;
  t475 = t445 + t455 + t462 + t466 + t470 + t474;
  t533 = -0.0265*t327;
  t534 = 0.0695*t356;
  t535 = t533 + t534;
  t503 = t465*t407;
  t504 = t378*t423;
  t507 = 0.0225*t506;
  t508 = t378*t399;
  t509 = t471 + t508;
  t510 = -0.0265*t509;
  t513 = -1.*t367*t444;
  t514 = -1.*t267*t454;
  t515 = t415*t407;
  t516 = t461*t423;
  t517 = 0.0225*t485;
  t518 = t461*t399;
  t519 = t434 + t518;
  t520 = -0.0265*t519;
  t521 = t513 + t514 + t515 + t516 + t517 + t520;
  t524 = t378*t478;
  t525 = t465*t482;
  t526 = -0.0265*t506;
  t527 = -1.*t380*t465;
  t528 = t527 + t435;
  t529 = 0.0225*t528;
  t530 = t524 + t525 + t526 + t529;
  t479 = t461*t478;
  t483 = t415*t482;
  t486 = -0.0265*t485;
  t487 = -1.*t380*t415;
  t488 = t487 + t472;
  t489 = 0.0225*t488;
  t490 = t479 + t483 + t486 + t489;
  t547 = -1.*t380*t461;
  t548 = t547 + t505;
  t493 = t415*t478;
  t494 = t378*t482;
  t495 = -1.*t380*t378;
  t496 = t495 + t484;
  t497 = 0.0225*t496;
  t498 = -0.0265*t436;
  t499 = t493 + t494 + t497 + t498;
  t545 = t465*t478;
  t546 = t461*t482;
  t549 = 0.0225*t548;
  t550 = -0.0265*t473;
  t551 = t545 + t546 + t549 + t550;
  t557 = -0.0265*t380;
  t558 = 0.2375*t399;
  t559 = t557 + t558;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=(-1.*t267*t363 - 1.*t367*t371 + t409 + t426 + t433 + t438)*var1[9] + t475*var1[10] + t490*var1[11];
  p_output1[5]=(t363*t367 - 1.*t267*t371 + t503 + t504 + t507 + t510)*var1[9] + t521*var1[10] + t499*var1[11];
  p_output1[6]=t475*var1[9] + (t409 + t426 + t433 + t438 + t367*t454 + t267*t535)*var1[10] + t530*var1[11];
  p_output1[7]=t521*var1[9] + (t267*t454 + t503 + t504 + t507 + t510 - 1.*t367*t535)*var1[10] + t551*var1[11];
  p_output1[8]=t490*var1[9] + t530*var1[10] + (t525 + 0.0225*(t518 + t527) - 0.0265*t548 + t461*t559)*var1[11];
  p_output1[9]=t499*var1[9] + t551*var1[10] + (-0.0265*t488 + t546 + 0.0225*(t429 + t547) + t415*t559)*var1[11];
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
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
    ( !(mrows == 14 && ncols == 1) && 
      !(mrows == 1 && ncols == 14))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "dJ_rightToe.hh"

namespace SymFunction
{

void dJ_rightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
