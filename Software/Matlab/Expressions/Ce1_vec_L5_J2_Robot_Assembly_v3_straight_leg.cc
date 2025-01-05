/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:30 GMT-05:00
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
  double t519;
  double t505;
  double t516;
  double t521;
  double t495;
  double t517;
  double t525;
  double t526;
  double t533;
  double t537;
  double t538;
  double t539;
  double t545;
  double t550;
  double t551;
  double t553;
  double t560;
  double t543;
  double t547;
  double t554;
  double t555;
  double t559;
  double t561;
  double t566;
  double t567;
  double t568;
  double t527;
  double t574;
  double t575;
  double t577;
  double t572;
  double t573;
  double t578;
  double t579;
  double t580;
  double t583;
  double t584;
  double t587;
  double t588;
  double t540;
  double t542;
  double t564;
  double t565;
  double t600;
  double t602;
  double t604;
  double t615;
  double t617;
  double t611;
  double t612;
  double t618;
  double t619;
  double t623;
  double t624;
  double t626;
  double t627;
  double t628;
  double t630;
  double t635;
  double t636;
  double t638;
  double t639;
  double t640;
  double t613;
  double t614;
  double t625;
  double t632;
  double t633;
  double t648;
  double t649;
  double t650;
  double t634;
  double t642;
  double t644;
  double t659;
  double t661;
  double t662;
  double t663;
  double t664;
  double t656;
  double t657;
  double t654;
  double t655;
  double t670;
  double t684;
  double t685;
  double t687;
  double t688;
  double t689;
  double t675;
  double t677;
  double t693;
  double t694;
  double t695;
  double t686;
  double t690;
  double t691;
  double t679;
  double t714;
  double t715;
  double t716;
  double t717;
  double t718;
  double t708;
  double t709;
  double t710;
  double t711;
  double t712;
  double t722;
  double t723;
  double t746;
  double t747;
  double t748;
  t519 = Cos(var1[5]);
  t505 = Cos(var1[6]);
  t516 = Sin(var1[5]);
  t521 = Sin(var1[6]);
  t495 = Sin(var1[2]);
  t517 = -1.*t505*t516;
  t525 = -1.*t519*t521;
  t526 = t517 + t525;
  t533 = Cos(var1[2]);
  t537 = t519*t505;
  t538 = -1.*t516*t521;
  t539 = t537 + t538;
  t545 = t533*t539;
  t550 = t505*t516;
  t551 = t519*t521;
  t553 = t550 + t551;
  t560 = -1.*t495*t539;
  t543 = t495*t526;
  t547 = t543 + t545;
  t554 = -1.*t495*t553;
  t555 = t554 + t545;
  t559 = t533*t526;
  t561 = t559 + t560;
  t566 = t533*t553;
  t567 = t495*t539;
  t568 = t566 + t567;
  t527 = -1.*t495*t526;
  t574 = -1.*t519*t505;
  t575 = t516*t521;
  t577 = t574 + t575;
  t572 = 0.19964*t547*t555;
  t573 = 0.19964*t561*t568;
  t578 = t533*t577;
  t579 = t527 + t578;
  t580 = 0.19964*t547*t579;
  t583 = t495*t577;
  t584 = t559 + t583;
  t587 = 0.19964*t561*t584;
  t588 = t572 + t573 + t580 + t587;
  t540 = -1.*t533*t539;
  t542 = t527 + t540;
  t564 = -1.*t533*t553;
  t565 = t564 + t560;
  t600 = 0.39928*t555*t561;
  t602 = 0.39928*t561*t579;
  t604 = t600 + t602;
  t615 = -1.*t505;
  t617 = 1. + t615;
  t611 = -1.*t519;
  t612 = 1. + t611;
  t618 = -0.0265*t617;
  t619 = -0.025226*t505;
  t623 = -0.07700600000000002*t521;
  t624 = t618 + t619 + t623;
  t626 = -0.2375*t617;
  t627 = -0.314506*t505;
  t628 = -0.0012740000000000008*t521;
  t630 = t626 + t627 + t628;
  t635 = -0.0695*t612;
  t636 = -0.0265*t516;
  t638 = -1.*t516*t624;
  t639 = t519*t630;
  t640 = t635 + t636 + t638 + t639;
  t613 = -0.0265*t612;
  t614 = 0.0695*t516;
  t625 = t519*t624;
  t632 = t516*t630;
  t633 = t613 + t614 + t625 + t632;
  t648 = -1.*t640*t526;
  t649 = -1.*t633*t539;
  t650 = t648 + t649;
  t634 = t633*t553;
  t642 = t640*t539;
  t644 = t634 + t642;
  t659 = -0.0265*t519;
  t661 = -0.0695*t516;
  t662 = -1.*t519*t624;
  t663 = -1.*t516*t630;
  t664 = t659 + t661 + t662 + t663;
  t656 = 0.0695*t519;
  t657 = t656 + t636 + t638 + t639;
  t654 = 0.19964*t561*t650;
  t655 = t640*t526;
  t670 = t633*t539;
  t684 = -0.07700600000000002*t505;
  t685 = t684 + t628;
  t687 = -0.0012740000000000008*t505;
  t688 = 0.07700600000000002*t521;
  t689 = t687 + t688;
  t675 = 0.19964*t644*t579;
  t677 = -1.*t633*t526;
  t693 = -1.*t516*t685;
  t694 = t519*t689;
  t695 = t693 + t694;
  t686 = t519*t685;
  t690 = t516*t689;
  t691 = t686 + t690;
  t679 = -1.*t640*t577;
  t714 = -0.0265*t505;
  t715 = -1.*t505*t624;
  t716 = 0.0695*t521;
  t717 = t630*t521;
  t718 = t714 + t715 + t716 + t717;
  t708 = 0.0695*t505;
  t709 = t505*t630;
  t710 = 0.0265*t521;
  t711 = t624*t521;
  t712 = t708 + t709 + t710 + t711;
  t722 = 0.19964*t718*t561;
  t723 = 0.19964*t712*t579;
  t746 = -0.0002543413600000002*t561;
  t747 = -0.015373477840000005*t579;
  t748 = t746 + t747;
  p_output1[0]=var2[1]*(-0.5*(0.19964*t542*t547 + 0.19964*Power(t555,2) + 0.19964*Power(t561,2) + 0.19964*t565*t568)*var2[2] - 0.5*t588*var2[5] - 0.5*t588*var2[6]);
  p_output1[1]=var2[1]*(-0.5*(0.39928*t542*t561 + 0.39928*t555*t565)*var2[2] - 0.5*t604*var2[5] - 0.5*t604*var2[6]);
  p_output1[2]=var2[1]*(-0.5*(0.19964*t542*t644 + 0.19964*t565*t650)*var2[2] - 0.5*(t654 + 0.19964*t561*(t655 + t553*t657 + t539*t664 + t670) + t675 + 0.19964*t555*(-1.*t539*t657 - 1.*t526*t664 + t677 + t679))*var2[5] - 0.5*(t654 + t675 + 0.19964*t555*(t677 + t679 - 1.*t539*t691 - 1.*t526*t695) + 0.19964*t561*(t655 + t670 + t553*t691 + t539*t695))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(-0.5*(0.19964*t542*t712 + 0.19964*t565*t718)*var2[2] - 0.5*(t722 + t723)*var2[5] - 0.5*(0.19964*t561*(0.0265*t505 - 0.0695*t521 + t505*t624 - 1.*t521*t630 + t521*t685 + t505*t689) + 0.19964*t555*(-1.*t505*t685 + t521*t689 + t708 + t709 + t710 + t711) + t722 + t723)*var2[6]);
  p_output1[6]=var2[1]*(-0.5*(-0.015373477840000005*t542 - 0.0002543413600000002*t565)*var2[2] - 0.5*t748*var2[5] - 0.5*t748*var2[6]);
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

#include "Ce1_vec_L5_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L5_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
