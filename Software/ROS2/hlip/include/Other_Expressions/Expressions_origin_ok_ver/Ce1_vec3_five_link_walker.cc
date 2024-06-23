/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:30 GMT-05:00
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
  double t306;
  double t372;
  double t307;
  double t369;
  double t471;
  double t486;
  double t503;
  double t504;
  double t505;
  double t472;
  double t487;
  double t490;
  double t510;
  double t511;
  double t514;
  double t515;
  double t516;
  double t502;
  double t506;
  double t507;
  double t517;
  double t520;
  double t521;
  double t523;
  double t526;
  double t527;
  double t568;
  double t570;
  double t571;
  double t531;
  double t556;
  double t559;
  double t583;
  double t584;
  double t585;
  double t522;
  double t528;
  double t529;
  double t530;
  double t561;
  double t566;
  double t579;
  double t580;
  double t581;
  double t582;
  double t586;
  double t587;
  double t588;
  double t589;
  double t590;
  double t611;
  double t631;
  double t633;
  double t616;
  double t635;
  double t636;
  double t622;
  double t366;
  double t373;
  double t374;
  double t375;
  double t376;
  double t379;
  double t460;
  double t466;
  double t467;
  double t667;
  double t669;
  double t678;
  double t680;
  double t684;
  double t685;
  double t691;
  double t697;
  double t700;
  double t701;
  double t702;
  double t703;
  double t708;
  double t709;
  double t710;
  double t704;
  double t705;
  double t706;
  double t679;
  double t681;
  double t682;
  double t711;
  double t712;
  double t713;
  double t668;
  double t670;
  double t671;
  double t672;
  double t673;
  double t674;
  double t675;
  double t676;
  double t677;
  double t695;
  double t707;
  double t714;
  double t715;
  double t726;
  double t727;
  double t720;
  double t721;
  double t722;
  double t717;
  double t729;
  double t730;
  double t731;
  double t738;
  double t739;
  double t740;
  double t728;
  double t733;
  double t734;
  double t741;
  double t742;
  double t743;
  double t744;
  double t745;
  double t746;
  double t755;
  double t756;
  double t748;
  double t758;
  double t759;
  double t750;
  double t662;
  double t663;
  double t560;
  double t572;
  double t573;
  double t651;
  double t652;
  double t615;
  double t617;
  double t625;
  double t775;
  double t634;
  double t639;
  double t642;
  double t777;
  double t778;
  double t779;
  double t780;
  double t644;
  double t645;
  double t646;
  double t771;
  double t772;
  double t773;
  double t774;
  double t718;
  double t795;
  double t796;
  double t797;
  double t798;
  double t719;
  double t732;
  double t735;
  double t736;
  double t802;
  double t683;
  double t696;
  double t747;
  double t749;
  double t751;
  double t808;
  double t757;
  double t760;
  double t761;
  double t810;
  double t811;
  double t812;
  double t763;
  double t764;
  double t765;
  double t840;
  double t841;
  double t842;
  double t843;
  double t845;
  double t846;
  double t847;
  double t869;
  double t870;
  double t871;
  double t872;
  double t874;
  double t875;
  double t876;
  t306 = Cos(var1[3]);
  t372 = Sin(var1[3]);
  t307 = Sin(var1[2]);
  t369 = Cos(var1[2]);
  t471 = Cos(var1[4]);
  t486 = Sin(var1[4]);
  t503 = t306*t471;
  t504 = -1.*t372*t486;
  t505 = t503 + t504;
  t472 = -1.*t471*t372;
  t487 = -1.*t306*t486;
  t490 = t472 + t487;
  t510 = -1.*t471;
  t511 = 1. + t510;
  t514 = 0.5*t511;
  t515 = 0.671885*t471;
  t516 = t514 + t515;
  t502 = t307*t490;
  t506 = t369*t505;
  t507 = t502 + t506;
  t517 = t516*t372;
  t520 = 0.171885*t306*t486;
  t521 = t517 + t520;
  t523 = t306*t516;
  t526 = -0.171885*t372*t486;
  t527 = t523 + t526;
  t568 = t471*t372;
  t570 = t306*t486;
  t571 = t568 + t570;
  t531 = -1.*t516*t372;
  t556 = -0.171885*t306*t486;
  t559 = t531 + t556;
  t583 = -1.*t306*t471;
  t584 = t372*t486;
  t585 = t583 + t584;
  t522 = -1.*t521*t505;
  t528 = -1.*t490*t527;
  t529 = t522 + t528;
  t530 = 0.85216*t507*t529;
  t561 = t521*t505;
  t566 = t490*t527;
  t579 = t521*t571;
  t580 = t505*t527;
  t581 = t579 + t580;
  t582 = t369*t490;
  t586 = t307*t585;
  t587 = t582 + t586;
  t588 = 0.85216*t581*t587;
  t589 = t369*t571;
  t590 = t307*t505;
  t611 = t589 + t590;
  t631 = -0.171885*t471*t372;
  t633 = t631 + t556;
  t616 = -1.*t490*t521;
  t635 = 0.171885*t306*t471;
  t636 = t635 + t526;
  t622 = -1.*t527*t585;
  t366 = -1.*t306*t307;
  t373 = -1.*t369*t372;
  t374 = t366 + t373;
  t375 = Power(t306,2);
  t376 = 0.1494*t375;
  t379 = Power(t372,2);
  t460 = 0.1494*t379;
  t466 = t376 + t460;
  t467 = 3.4261*t374*t466;
  t667 = Cos(var1[5]);
  t669 = Sin(var1[5]);
  t678 = Cos(var1[6]);
  t680 = Sin(var1[6]);
  t684 = t667*t678;
  t685 = -1.*t669*t680;
  t691 = t684 + t685;
  t697 = -1.*t678;
  t700 = 1. + t697;
  t701 = 0.5*t700;
  t702 = 0.671885*t678;
  t703 = t701 + t702;
  t708 = -1.*t678*t669;
  t709 = -1.*t667*t680;
  t710 = t708 + t709;
  t704 = t703*t669;
  t705 = 0.171885*t667*t680;
  t706 = t704 + t705;
  t679 = t678*t669;
  t681 = t667*t680;
  t682 = t679 + t681;
  t711 = t667*t703;
  t712 = -0.171885*t669*t680;
  t713 = t711 + t712;
  t668 = -1.*t667*t307;
  t670 = -1.*t369*t669;
  t671 = t668 + t670;
  t672 = Power(t667,2);
  t673 = 0.1494*t672;
  t674 = Power(t669,2);
  t675 = 0.1494*t674;
  t676 = t673 + t675;
  t677 = 3.4261*t671*t676;
  t695 = t369*t691;
  t707 = -1.*t706*t691;
  t714 = -1.*t710*t713;
  t715 = t707 + t714;
  t726 = t307*t710;
  t727 = t726 + t695;
  t720 = t706*t682;
  t721 = t691*t713;
  t722 = t720 + t721;
  t717 = t369*t710;
  t729 = -1.*t703*t669;
  t730 = -0.171885*t667*t680;
  t731 = t729 + t730;
  t738 = -1.*t667*t678;
  t739 = t669*t680;
  t740 = t738 + t739;
  t728 = 0.85216*t727*t715;
  t733 = t706*t691;
  t734 = t710*t713;
  t741 = t307*t740;
  t742 = t717 + t741;
  t743 = 0.85216*t722*t742;
  t744 = t369*t682;
  t745 = t307*t691;
  t746 = t744 + t745;
  t755 = -0.171885*t678*t669;
  t756 = t755 + t730;
  t748 = -1.*t710*t706;
  t758 = 0.171885*t667*t678;
  t759 = t758 + t712;
  t750 = -1.*t713*t740;
  t662 = -1.*t307*t505;
  t663 = t582 + t662;
  t560 = t559*t505;
  t572 = t571*t527;
  t573 = t560 + t561 + t566 + t572;
  t651 = -1.*t307*t571;
  t652 = t651 + t506;
  t615 = -1.*t490*t559;
  t617 = -1.*t505*t527;
  t625 = t615 + t616 + t617 + t622;
  t775 = 0.85216*t663*t529;
  t634 = t633*t505;
  t639 = t571*t636;
  t642 = t634 + t561 + t566 + t639;
  t777 = -1.*t307*t490;
  t778 = t369*t585;
  t779 = t777 + t778;
  t780 = 0.85216*t581*t779;
  t644 = -1.*t490*t633;
  t645 = -1.*t505*t636;
  t646 = t644 + t616 + t645 + t622;
  t771 = -1.*t369*t306;
  t772 = t307*t372;
  t773 = t771 + t772;
  t774 = 3.4261*t773*t466;
  t718 = -1.*t307*t691;
  t795 = -1.*t369*t667;
  t796 = t307*t669;
  t797 = t795 + t796;
  t798 = 3.4261*t797*t676;
  t719 = t717 + t718;
  t732 = t731*t691;
  t735 = t682*t713;
  t736 = t732 + t733 + t734 + t735;
  t802 = -1.*t307*t710;
  t683 = -1.*t307*t682;
  t696 = t683 + t695;
  t747 = -1.*t710*t731;
  t749 = -1.*t691*t713;
  t751 = t747 + t748 + t749 + t750;
  t808 = 0.85216*t719*t715;
  t757 = t756*t691;
  t760 = t682*t759;
  t761 = t757 + t733 + t734 + t760;
  t810 = t369*t740;
  t811 = t802 + t810;
  t812 = 0.85216*t722*t811;
  t763 = -1.*t710*t756;
  t764 = -1.*t691*t759;
  t765 = t763 + t748 + t764 + t750;
  t840 = t516*t471;
  t841 = Power(t486,2);
  t842 = 0.171885*t841;
  t843 = t840 + t842;
  t845 = t516*t486;
  t846 = -0.171885*t471*t486;
  t847 = t845 + t846;
  t869 = t703*t678;
  t870 = Power(t680,2);
  t871 = 0.171885*t870;
  t872 = t869 + t871;
  t874 = t703*t680;
  t875 = -0.171885*t678*t680;
  t876 = t874 + t875;
  p_output1[0]=var2[2]*(-0.5*(-3.70591*t307 + t467 + 0.85216*t529*t652 + 0.85216*t581*t663 + t677 + 0.85216*t696*t715 + 0.85216*t719*t722)*var2[2] - 0.5*(t467 + t530 + 0.85216*t507*t573 + t588 + 0.85216*t611*t625)*var2[3] - 0.5*(t530 + t588 + 0.85216*t507*t642 + 0.85216*t611*t646)*var2[4] - 0.5*(t677 + t728 + 0.85216*t727*t736 + t743 + 0.85216*t746*t751)*var2[5] - 0.5*(t728 + t743 + 0.85216*t727*t761 + 0.85216*t746*t765)*var2[6]);
  p_output1[1]=var2[2]*(-0.5*(-3.70591*t369 + 0.85216*t529*(-1.*t369*t571 + t662) + 0.85216*t715*(-1.*t369*t682 + t718) + t774 + 0.85216*t581*(-1.*t369*t505 + t777) + t798 + 0.85216*t722*(-1.*t369*t691 + t802))*var2[2] - 0.5*(0.85216*t625*t652 + 0.85216*t573*t663 + t774 + t775 + t780)*var2[3] - 0.5*(0.85216*t646*t652 + 0.85216*t642*t663 + t775 + t780)*var2[4] - 0.5*(0.85216*t719*t736 + 0.85216*t696*t751 + t798 + t808 + t812)*var2[5] - 0.5*(0.85216*t719*t761 + 0.85216*t696*t765 + t808 + t812)*var2[6]);
  p_output1[2]=var2[2]*(-0.5*(1.70432*t573*t581 + 1.70432*t529*t625)*var2[3] - 0.5*(1.70432*t581*t642 + 1.70432*t529*t646)*var2[4] - 0.5*(1.70432*t722*t736 + 1.70432*t715*t751)*var2[5] - 0.5*(1.70432*t722*t761 + 1.70432*t715*t765)*var2[6]);
  p_output1[3]=var2[2]*(-0.5*(0.85216*t573*t843 + 0.85216*t625*t847)*var2[3] - 0.5*(0.85216*(0.171885*t471*t486 - 1.*t486*t516)*t581 + 0.85216*t529*(-0.171885*Power(t471,2) + t840) + 0.85216*t642*t843 + 0.85216*t646*t847)*var2[4]);
  p_output1[4]=var2[2]*(-0.0732367608*t573*var2[3] - 0.0732367608*t642*var2[4]);
  p_output1[5]=var2[2]*(-0.5*(0.85216*t736*t872 + 0.85216*t751*t876)*var2[5] - 0.5*(0.85216*(0.171885*t678*t680 - 1.*t680*t703)*t722 + 0.85216*t715*(-0.171885*Power(t678,2) + t869) + 0.85216*t761*t872 + 0.85216*t765*t876)*var2[6]);
  p_output1[6]=var2[2]*(-0.0732367608*t736*var2[5] - 0.0732367608*t761*var2[6]);
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

#include "Ce1_vec3_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
