/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:38:09 GMT-05:00
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
  double t453;
  double t475;
  double t522;
  double t518;
  double t486;
  double t487;
  double t488;
  double t499;
  double t503;
  double t504;
  double t460;
  double t478;
  double t479;
  double t531;
  double t516;
  double t535;
  double t536;
  double t538;
  double t482;
  double t507;
  double t512;
  double t563;
  double t565;
  double t567;
  double t542;
  double t556;
  double t557;
  double t560;
  double t568;
  double t588;
  double t589;
  double t590;
  double t621;
  double t622;
  double t611;
  double t614;
  double t617;
  double t520;
  double t525;
  double t529;
  double t591;
  double t592;
  double t595;
  double t597;
  double t600;
  double t604;
  double t569;
  double t574;
  double t640;
  double t629;
  double t630;
  double t631;
  double t530;
  double t544;
  double t650;
  double t610;
  double t618;
  double t620;
  double t624;
  double t627;
  double t651;
  double t652;
  double t653;
  double t665;
  double t666;
  double t667;
  double t670;
  double t671;
  double t672;
  double t674;
  double t675;
  double t676;
  double t668;
  double t677;
  double t681;
  double t700;
  double t701;
  double t702;
  double t696;
  double t697;
  double t698;
  double t683;
  t453 = Cos(var1[6]);
  t475 = Sin(var1[6]);
  t522 = Cos(var1[5]);
  t518 = Sin(var1[5]);
  t486 = -1.*t453;
  t487 = 1. + t486;
  t488 = -0.16*t487;
  t499 = -0.167368*t453;
  t503 = 0.022659*t475;
  t504 = t488 + t499 + t503;
  t460 = -0.022659*t453;
  t478 = -0.007367999999999986*t475;
  t479 = t460 + t478;
  t531 = Cos(var1[2]);
  t516 = Sin(var1[2]);
  t535 = t522*t453;
  t536 = -1.*t518*t475;
  t538 = t535 + t536;
  t482 = -1.*t453*t479;
  t507 = t504*t475;
  t512 = t482 + t507;
  t563 = -1.*t453*t518;
  t565 = -1.*t522*t475;
  t567 = t563 + t565;
  t542 = t531*t538;
  t556 = t453*t504;
  t557 = t479*t475;
  t560 = t556 + t557;
  t568 = t531*t567;
  t588 = t516*t567;
  t589 = t588 + t542;
  t590 = 0.14994*t512*t589;
  t621 = -0.007367999999999986*t453;
  t622 = t621 + t503;
  t611 = 0.022659*t453;
  t614 = 0.007367999999999986*t475;
  t617 = t611 + t614;
  t520 = t453*t518;
  t525 = t522*t475;
  t529 = t520 + t525;
  t591 = -1.*t522*t453;
  t592 = t518*t475;
  t595 = t591 + t592;
  t597 = t516*t595;
  t600 = t568 + t597;
  t604 = 0.14994*t560*t600;
  t569 = -1.*t516*t538;
  t574 = t568 + t569;
  t640 = -1.*t516*t567;
  t629 = -1.*t453*t622;
  t630 = t617*t475;
  t631 = t556 + t629 + t557 + t630;
  t530 = -1.*t516*t529;
  t544 = t530 + t542;
  t650 = 0.14994*t512*t574;
  t610 = t453*t479;
  t618 = t453*t617;
  t620 = -1.*t504*t475;
  t624 = t622*t475;
  t627 = t610 + t618 + t620 + t624;
  t651 = t531*t595;
  t652 = t640 + t651;
  t653 = 0.14994*t560*t652;
  t665 = -1.*t518*t479;
  t666 = t522*t504;
  t667 = t665 + t666;
  t670 = -1.*t522*t479;
  t671 = -1.*t518*t504;
  t672 = t670 + t671;
  t674 = t522*t479;
  t675 = t518*t504;
  t676 = t674 + t675;
  t668 = t667*t567;
  t677 = t676*t538;
  t681 = -1.*t676*t567;
  t700 = t522*t617;
  t701 = -1.*t518*t622;
  t702 = t700 + t701;
  t696 = t518*t617;
  t697 = t522*t622;
  t698 = t696 + t697;
  t683 = -1.*t667*t595;
  p_output1[0]=var2[5]*(-0.5*(0.14994*t512*t544 + 0.14994*t560*t574)*var2[2] - 0.5*(t590 + t604)*var2[5] - 0.5*(t590 + t604 + 0.14994*t589*t627 + 0.14994*(t529*t531 + t516*t538)*t631)*var2[6]);
  p_output1[1]=var2[5]*(-0.5*(0.14994*t512*(-1.*t529*t531 + t569) + 0.14994*t560*(-1.*t531*t538 + t640))*var2[2] - 0.5*(t650 + t653)*var2[5] - 0.5*(0.14994*t574*t627 + 0.14994*t544*t631 + t650 + t653)*var2[6]);
  p_output1[2]=var2[5]*(-0.5*(0.14994*t560*(t529*t667 + t668 + t538*t672 + t677) + 0.14994*t512*(-1.*t538*t667 - 1.*t567*t672 + t681 + t683))*var2[5] - 0.5*(0.14994*t627*(t538*t667 + t529*t676) + 0.14994*t631*(-1.*t567*t667 - 1.*t538*t676) + 0.14994*t560*(t668 + t677 + t529*t698 + t538*t702) + 0.14994*t512*(t681 + t683 - 1.*t538*t698 - 1.*t567*t702))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(0.29988*t560*t627 + 0.29988*t512*t631)*var2[5]*var2[6];
  p_output1[6]=-0.5*(-0.0011047579199999977*t627 + 0.0033974904599999994*t631)*var2[5]*var2[6];
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

#include "Ce1_vec_L5_J6_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L5_J6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
