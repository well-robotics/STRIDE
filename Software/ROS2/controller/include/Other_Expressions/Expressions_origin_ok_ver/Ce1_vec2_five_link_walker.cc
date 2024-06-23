/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:29 GMT-05:00
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
  double t154;
  double t74;
  double t142;
  double t166;
  double t73;
  double t242;
  double t252;
  double t253;
  double t254;
  double t255;
  double t148;
  double t217;
  double t240;
  double t257;
  double t265;
  double t266;
  double t241;
  double t256;
  double t279;
  double t280;
  double t281;
  double t287;
  double t298;
  double t299;
  double t308;
  double t309;
  double t316;
  double t349;
  double t350;
  double t351;
  double t276;
  double t277;
  double t278;
  double t282;
  double t283;
  double t284;
  double t285;
  double t286;
  double t300;
  double t301;
  double t302;
  double t303;
  double t304;
  double t305;
  double t317;
  double t318;
  double t319;
  double t320;
  double t321;
  double t346;
  double t356;
  double t358;
  double t360;
  double t361;
  double t362;
  double t363;
  double t380;
  double t401;
  double t405;
  double t406;
  double t407;
  double t423;
  double t424;
  double t425;
  double t436;
  double t438;
  double t437;
  double t439;
  double t440;
  double t442;
  double t447;
  double t452;
  double t458;
  double t461;
  double t462;
  double t463;
  double t469;
  double t457;
  double t459;
  double t464;
  double t465;
  double t468;
  double t470;
  double t475;
  double t480;
  double t485;
  double t441;
  double t493;
  double t494;
  double t495;
  double t412;
  double t415;
  double t416;
  double t417;
  double t419;
  double t421;
  double t426;
  double t429;
  double t432;
  double t433;
  double t434;
  double t435;
  double t491;
  double t492;
  double t496;
  double t497;
  double t498;
  double t499;
  double t500;
  double t501;
  double t508;
  double t509;
  double t512;
  double t513;
  double t370;
  double t371;
  double t377;
  double t378;
  double t453;
  double t454;
  double t473;
  double t474;
  double t518;
  double t519;
  double t524;
  double t525;
  double t538;
  double t539;
  double t540;
  double t541;
  double t542;
  double t543;
  double t544;
  double t545;
  double t547;
  double t548;
  double t549;
  double t553;
  double t554;
  double t555;
  double t546;
  double t550;
  double t551;
  double t552;
  double t557;
  double t558;
  double t562;
  double t563;
  double t564;
  double t565;
  double t574;
  double t575;
  double t567;
  double t577;
  double t578;
  double t569;
  double t532;
  double t533;
  double t534;
  double t535;
  double t536;
  double t537;
  double t597;
  double t598;
  double t599;
  double t600;
  double t601;
  double t602;
  double t603;
  double t604;
  double t606;
  double t607;
  double t608;
  double t591;
  double t592;
  double t593;
  double t594;
  double t595;
  double t596;
  double t605;
  double t609;
  double t610;
  double t612;
  double t613;
  double t614;
  double t619;
  double t620;
  double t621;
  double t618;
  double t623;
  double t624;
  double t628;
  double t637;
  double t638;
  double t630;
  double t640;
  double t641;
  double t632;
  double t653;
  double t659;
  double t660;
  double t661;
  double t654;
  double t655;
  double t656;
  double t657;
  double t665;
  double t666;
  double t686;
  double t692;
  double t693;
  double t694;
  double t687;
  double t688;
  double t689;
  double t690;
  double t698;
  double t699;
  t154 = Cos(var1[3]);
  t74 = Cos(var1[4]);
  t142 = Sin(var1[3]);
  t166 = Sin(var1[4]);
  t73 = Sin(var1[2]);
  t242 = Cos(var1[2]);
  t252 = t154*t74;
  t253 = -1.*t142*t166;
  t254 = t252 + t253;
  t255 = t242*t254;
  t148 = -1.*t74*t142;
  t217 = -1.*t154*t166;
  t240 = t148 + t217;
  t257 = t74*t142;
  t265 = t154*t166;
  t266 = t257 + t265;
  t241 = t73*t240;
  t256 = t241 + t255;
  t279 = t242*t240;
  t280 = -1.*t73*t254;
  t281 = t279 + t280;
  t287 = -1.*t154*t74;
  t298 = t142*t166;
  t299 = t287 + t298;
  t308 = -1.*t154*t73;
  t309 = -1.*t242*t142;
  t316 = t308 + t309;
  t349 = t242*t154;
  t350 = -1.*t73*t142;
  t351 = t349 + t350;
  t276 = -1.*t73*t266;
  t277 = t276 + t255;
  t278 = 0.85216*t256*t277;
  t282 = t242*t266;
  t283 = t73*t254;
  t284 = t282 + t283;
  t285 = 0.85216*t281*t284;
  t286 = -1.*t73*t240;
  t300 = t242*t299;
  t301 = t286 + t300;
  t302 = 0.85216*t256*t301;
  t303 = t73*t299;
  t304 = t279 + t303;
  t305 = 0.85216*t281*t304;
  t317 = Power(t316,2);
  t318 = 3.4261*t317;
  t319 = t154*t73;
  t320 = t242*t142;
  t321 = t319 + t320;
  t346 = 3.4261*t316*t321;
  t356 = Power(t351,2);
  t358 = 3.4261*t356;
  t360 = -1.*t242*t154;
  t361 = t73*t142;
  t362 = t360 + t361;
  t363 = 3.4261*t351*t362;
  t380 = Cos(var1[5]);
  t401 = -1.*t380*t73;
  t405 = Sin(var1[5]);
  t406 = -1.*t242*t405;
  t407 = t401 + t406;
  t423 = t242*t380;
  t424 = -1.*t73*t405;
  t425 = t423 + t424;
  t436 = Cos(var1[6]);
  t438 = Sin(var1[6]);
  t437 = -1.*t436*t405;
  t439 = -1.*t380*t438;
  t440 = t437 + t439;
  t442 = t380*t436;
  t447 = -1.*t405*t438;
  t452 = t442 + t447;
  t458 = t242*t452;
  t461 = t436*t405;
  t462 = t380*t438;
  t463 = t461 + t462;
  t469 = -1.*t73*t452;
  t457 = t73*t440;
  t459 = t457 + t458;
  t464 = -1.*t73*t463;
  t465 = t464 + t458;
  t468 = t242*t440;
  t470 = t468 + t469;
  t475 = t242*t463;
  t480 = t73*t452;
  t485 = t475 + t480;
  t441 = -1.*t73*t440;
  t493 = -1.*t380*t436;
  t494 = t405*t438;
  t495 = t493 + t494;
  t412 = Power(t407,2);
  t415 = 3.4261*t412;
  t416 = t380*t73;
  t417 = t242*t405;
  t419 = t416 + t417;
  t421 = 3.4261*t407*t419;
  t426 = Power(t425,2);
  t429 = 3.4261*t426;
  t432 = -1.*t242*t380;
  t433 = t73*t405;
  t434 = t432 + t433;
  t435 = 3.4261*t425*t434;
  t491 = 0.85216*t459*t465;
  t492 = 0.85216*t470*t485;
  t496 = t242*t495;
  t497 = t441 + t496;
  t498 = 0.85216*t459*t497;
  t499 = t73*t495;
  t500 = t468 + t499;
  t501 = 0.85216*t470*t500;
  t508 = 1.70432*t277*t281;
  t509 = 1.70432*t281*t301;
  t512 = 6.8522*t316*t351;
  t513 = 6.8522*t316*t362;
  t370 = -1.*t242*t254;
  t371 = t286 + t370;
  t377 = -1.*t242*t266;
  t378 = t377 + t280;
  t453 = -1.*t242*t452;
  t454 = t441 + t453;
  t473 = -1.*t242*t463;
  t474 = t473 + t469;
  t518 = 6.8522*t407*t425;
  t519 = 6.8522*t407*t434;
  t524 = 1.70432*t465*t470;
  t525 = 1.70432*t470*t497;
  t538 = -1.*t74;
  t539 = 1. + t538;
  t540 = 0.5*t539;
  t541 = 0.671885*t74;
  t542 = t540 + t541;
  t543 = t542*t142;
  t544 = 0.171885*t154*t166;
  t545 = t543 + t544;
  t547 = t154*t542;
  t548 = -0.171885*t142*t166;
  t549 = t547 + t548;
  t553 = -1.*t542*t142;
  t554 = -0.171885*t154*t166;
  t555 = t553 + t554;
  t546 = -1.*t545*t254;
  t550 = -1.*t240*t549;
  t551 = t546 + t550;
  t552 = 0.85216*t281*t551;
  t557 = t545*t254;
  t558 = t240*t549;
  t562 = t545*t266;
  t563 = t254*t549;
  t564 = t562 + t563;
  t565 = 0.85216*t564*t301;
  t574 = -0.171885*t74*t142;
  t575 = t574 + t554;
  t567 = -1.*t240*t545;
  t577 = 0.171885*t154*t74;
  t578 = t577 + t548;
  t569 = -1.*t549*t299;
  t532 = Power(t154,2);
  t533 = 0.1494*t532;
  t534 = Power(t142,2);
  t535 = 0.1494*t534;
  t536 = t533 + t535;
  t537 = 3.4261*t362*t536;
  t597 = -1.*t436;
  t598 = 1. + t597;
  t599 = 0.5*t598;
  t600 = 0.671885*t436;
  t601 = t599 + t600;
  t602 = t601*t405;
  t603 = 0.171885*t380*t438;
  t604 = t602 + t603;
  t606 = t380*t601;
  t607 = -0.171885*t405*t438;
  t608 = t606 + t607;
  t591 = Power(t380,2);
  t592 = 0.1494*t591;
  t593 = Power(t405,2);
  t594 = 0.1494*t593;
  t595 = t592 + t594;
  t596 = 3.4261*t434*t595;
  t605 = -1.*t604*t452;
  t609 = -1.*t440*t608;
  t610 = t605 + t609;
  t612 = t604*t463;
  t613 = t452*t608;
  t614 = t612 + t613;
  t619 = -1.*t601*t405;
  t620 = -0.171885*t380*t438;
  t621 = t619 + t620;
  t618 = 0.85216*t470*t610;
  t623 = t604*t452;
  t624 = t440*t608;
  t628 = 0.85216*t614*t497;
  t637 = -0.171885*t436*t405;
  t638 = t637 + t620;
  t630 = -1.*t440*t604;
  t640 = 0.171885*t380*t436;
  t641 = t640 + t607;
  t632 = -1.*t608*t495;
  t653 = 0.51185934*t362;
  t659 = t542*t166;
  t660 = -0.171885*t74*t166;
  t661 = t659 + t660;
  t654 = t542*t74;
  t655 = Power(t166,2);
  t656 = 0.171885*t655;
  t657 = t654 + t656;
  t665 = 0.85216*t661*t281;
  t666 = 0.85216*t657*t301;
  t686 = 0.51185934*t434;
  t692 = t601*t438;
  t693 = -0.171885*t436*t438;
  t694 = t692 + t693;
  t687 = t601*t436;
  t688 = Power(t438,2);
  t689 = 0.171885*t688;
  t690 = t687 + t689;
  t698 = 0.85216*t694*t470;
  t699 = 0.85216*t690*t497;
  p_output1[0]=var2[1]*(-0.5*(0.85216*Power(t277,2) + 0.85216*Power(t281,2) + t318 + t346 + t358 + t363 + 0.85216*t256*t371 + 0.85216*t284*t378 + t415 + t421 + t429 + t435 + 0.85216*t454*t459 + 0.85216*Power(t465,2) + 0.85216*Power(t470,2) + 0.85216*t474*t485)*var2[2] - 0.5*(t278 + t285 + t302 + t305 + t318 + t346 + t358 + t363)*var2[3] - 0.5*(t278 + t285 + t302 + t305)*var2[4] - 0.5*(t415 + t421 + t429 + t435 + t491 + t492 + t498 + t501)*var2[5] - 0.5*(t491 + t492 + t498 + t501)*var2[6]);
  p_output1[1]=var2[1]*(-0.5*(1.70432*t281*t371 + 1.70432*t277*t378 + 1.70432*t454*t470 + 1.70432*t465*t474 + t512 + t513 + t518 + t519)*var2[2] - 0.5*(t508 + t509 + t512 + t513)*var2[3] - 0.5*(t508 + t509)*var2[4] - 0.5*(t518 + t519 + t524 + t525)*var2[5] - 0.5*(t524 + t525)*var2[6]);
  p_output1[2]=var2[1]*(-0.5*(-3.70591*t242 + t537 + 0.85216*t378*t551 + 0.85216*t371*t564 + t596 + 0.85216*t474*t610 + 0.85216*t454*t614)*var2[2] - 0.5*(t537 + t552 + 0.85216*t281*(t266*t549 + t254*t555 + t557 + t558) + t565 + 0.85216*t277*(-1.*t254*t549 - 1.*t240*t555 + t567 + t569))*var2[3] - 0.5*(t552 + t565 + 0.85216*t277*(t567 + t569 - 1.*t240*t575 - 1.*t254*t578) + 0.85216*t281*(t557 + t558 + t254*t575 + t266*t578))*var2[4] - 0.5*(t596 + t618 + 0.85216*t470*(t463*t608 + t452*t621 + t623 + t624) + t628 + 0.85216*t465*(-1.*t452*t608 - 1.*t440*t621 + t630 + t632))*var2[5] - 0.5*(t618 + t628 + 0.85216*t465*(t630 + t632 - 1.*t440*t638 - 1.*t452*t641) + 0.85216*t470*(t623 + t624 + t452*t638 + t463*t641))*var2[6]);
  p_output1[3]=var2[1]*(-0.5*(t653 + 0.85216*t371*t657 + 0.85216*t378*t661)*var2[2] - 0.5*(t653 + t665 + t666)*var2[3] - 0.5*(t665 + t666 + 0.85216*t281*(-1.*t166*t542 + 0.171885*t166*t74) + 0.85216*t277*(t654 - 0.171885*Power(t74,2)))*var2[4]);
  p_output1[4]=var2[1]*(-0.0732367608*t371*var2[2] - 0.0732367608*t301*var2[3] - 0.0732367608*t301*var2[4]);
  p_output1[5]=var2[1]*(-0.5*(t686 + 0.85216*t454*t690 + 0.85216*t474*t694)*var2[2] - 0.5*(t686 + t698 + t699)*var2[5] - 0.5*(0.85216*t470*(0.171885*t436*t438 - 1.*t438*t601) + 0.85216*t465*(-0.171885*Power(t436,2) + t687) + t698 + t699)*var2[6]);
  p_output1[6]=var2[1]*(-0.0732367608*t454*var2[2] - 0.0732367608*t497*var2[5] - 0.0732367608*t497*var2[6]);
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

#include "Ce1_vec2_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
