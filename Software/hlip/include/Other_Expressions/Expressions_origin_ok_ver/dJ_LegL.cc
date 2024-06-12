/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:45 GMT-05:00
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
  double t458;
  double t482;
  double t443;
  double t506;
  double t560;
  double t563;
  double t565;
  double t569;
  double t583;
  double t587;
  double t593;
  double t548;
  double t550;
  double t559;
  double t486;
  double t648;
  double t649;
  double t568;
  double t570;
  double t573;
  double t605;
  double t606;
  double t607;
  double t660;
  double t662;
  double t663;
  double t627;
  double t478;
  double t492;
  double t515;
  double t522;
  double t542;
  double t664;
  double t665;
  double t666;
  double t667;
  double t668;
  double t669;
  double t670;
  double t671;
  double t672;
  double t651;
  double t652;
  double t655;
  double t656;
  double t658;
  double t684;
  double t685;
  double t686;
  double t623;
  double t692;
  double t681;
  double t682;
  double t687;
  double t688;
  double t689;
  double t690;
  double t691;
  double t693;
  double t694;
  double t695;
  double t696;
  double t654;
  double t659;
  double t673;
  double t699;
  double t701;
  double t702;
  double t705;
  double t710;
  double t712;
  double t714;
  double t720;
  double t721;
  double t723;
  double t724;
  double t726;
  double t727;
  double t728;
  double t730;
  double t731;
  double t732;
  double t735;
  double t624;
  double t742;
  double t745;
  double t749;
  double t752;
  double t753;
  double t766;
  double t767;
  double t768;
  double t769;
  double t770;
  double t775;
  double t776;
  double t748;
  double t754;
  double t758;
  double t759;
  double t760;
  double t762;
  double t763;
  double t675;
  double t677;
  double t679;
  double t494;
  double t544;
  double t580;
  double t614;
  double t625;
  double t626;
  double t631;
  double t641;
  double t643;
  double t646;
  double t796;
  double t824;
  double t825;
  double t826;
  double t828;
  double t839;
  double t840;
  double t841;
  double t843;
  double t836;
  double t837;
  double t838;
  double t847;
  double t848;
  double t849;
  double t827;
  double t829;
  double t830;
  double t832;
  double t833;
  double t834;
  double t842;
  double t844;
  double t845;
  double t850;
  double t851;
  double t852;
  double t854;
  double t866;
  double t867;
  double t868;
  double t831;
  double t835;
  double t846;
  double t853;
  double t855;
  double t856;
  double t857;
  double t858;
  double t859;
  double t860;
  double t861;
  double t862;
  double t871;
  double t872;
  double t882;
  double t883;
  double t884;
  double t886;
  double t887;
  double t864;
  double t865;
  double t869;
  double t870;
  double t873;
  double t874;
  double t875;
  double t876;
  double t877;
  double t878;
  double t905;
  double t906;
  double t907;
  double t898;
  double t895;
  double t896;
  double t897;
  double t899;
  double t900;
  double t901;
  double t902;
  double t909;
  double t910;
  double t885;
  double t888;
  double t889;
  double t890;
  double t891;
  double t892;
  double t893;
  double t863;
  double t879;
  double t880;
  double t918;
  double t919;
  double t920;
  double t921;
  double t922;
  double t923;
  double t924;
  double t925;
  double t926;
  double t933;
  double t938;
  double t939;
  double t940;
  double t942;
  double t943;
  double t934;
  double t935;
  double t936;
  double t881;
  double t941;
  double t944;
  double t945;
  double t913;
  double t947;
  double t948;
  double t949;
  double t958;
  double t959;
  double t960;
  double t961;
  double t962;
  double t963;
  double t964;
  double t674;
  double t697;
  double t698;
  double t765;
  double t777;
  double t778;
  double t991;
  double t992;
  double t993;
  double t794;
  double t798;
  double t801;
  double t804;
  double t805;
  double t807;
  double t808;
  double t810;
  double t811;
  double t1002;
  double t1003;
  double t1004;
  double t1005;
  double t1006;
  double t1007;
  double t1008;
  double t1009;
  double t1010;
  double t1011;
  double t1012;
  double t781;
  double t782;
  double t783;
  double t784;
  double t785;
  double t787;
  double t788;
  double t789;
  double t791;
  double t792;
  double t795;
  double t1026;
  double t1027;
  double t1028;
  double t903;
  double t904;
  double t908;
  double t911;
  double t912;
  double t914;
  double t915;
  double t916;
  double t917;
  double t930;
  double t931;
  double t932;
  double t946;
  double t950;
  double t951;
  double t956;
  double t1053;
  double t1054;
  double t1055;
  double t955;
  double t957;
  double t965;
  double t966;
  double t967;
  double t968;
  double t969;
  double t1075;
  double t1076;
  double t1077;
  double t1044;
  double t1045;
  double t1046;
  double t894;
  double t927;
  double t928;
  t458 = Cos(var1[3]);
  t482 = Sin(var1[3]);
  t443 = Sin(var1[2]);
  t506 = Cos(var1[2]);
  t560 = Cos(var1[4]);
  t563 = -1.*t560;
  t565 = 1. + t563;
  t569 = Sin(var1[4]);
  t583 = -1.*t506*t458;
  t587 = -1.*t443*t482;
  t593 = t583 + t587;
  t548 = t458*t443;
  t550 = -1.*t506*t482;
  t559 = t548 + t550;
  t486 = -0.0695*t482;
  t648 = -1.*t458;
  t649 = 1. + t648;
  t568 = -0.0265*t565;
  t570 = -0.2375*t569;
  t573 = t568 + t570;
  t605 = -0.2375*t565;
  t606 = 0.0265*t569;
  t607 = t605 + t606;
  t660 = t506*t458;
  t662 = t443*t482;
  t663 = t660 + t662;
  t627 = t560*t559;
  t478 = 0.0265*t458;
  t492 = t478 + t486;
  t515 = -0.0695*t458;
  t522 = -0.0265*t482;
  t542 = t515 + t522;
  t664 = -1.*t663*t573;
  t665 = -1.*t559*t607;
  t666 = t560*t663;
  t667 = t559*t569;
  t668 = t666 + t667;
  t669 = 0.0265*t668;
  t670 = -1.*t663*t569;
  t671 = t627 + t670;
  t672 = -0.0225*t671;
  t651 = -0.0265*t649;
  t652 = t651 + t486;
  t655 = -0.0695*t649;
  t656 = 0.0265*t482;
  t658 = t655 + t656;
  t684 = -1.*t458*t443;
  t685 = t506*t482;
  t686 = t684 + t685;
  t623 = t560*t593;
  t692 = t560*t686;
  t681 = t443*t652;
  t682 = -1.*t506*t658;
  t687 = -1.*t686*t573;
  t688 = -1.*t663*t607;
  t689 = -1.*t686*t569;
  t690 = t666 + t689;
  t691 = -0.0225*t690;
  t693 = t663*t569;
  t694 = t692 + t693;
  t695 = 0.0265*t694;
  t696 = t681 + t682 + t687 + t688 + t691 + t695;
  t654 = -1.*t506*t652;
  t659 = -1.*t443*t658;
  t673 = t654 + t659 + t664 + t665 + t669 + t672;
  t699 = t506*t652;
  t701 = t443*t658;
  t702 = -1.*t593*t573;
  t705 = -1.*t686*t607;
  t710 = t686*t569;
  t712 = t623 + t710;
  t714 = 0.0265*t712;
  t720 = -1.*t593*t569;
  t721 = t692 + t720;
  t723 = -0.0225*t721;
  t724 = t699 + t701 + t702 + t705 + t714 + t723;
  t726 = 2.*t724*t696;
  t727 = 2.*t673*t696;
  t728 = t726 + t727;
  t730 = Power(t673,2);
  t731 = Power(t696,2);
  t732 = 0.00085849 + t730 + t731;
  t735 = Power(t732,-1.5);
  t624 = -1.*t559*t569;
  t742 = 0.0265*t560;
  t745 = t742 + t570;
  t749 = -0.2375*t560;
  t752 = -0.0265*t569;
  t753 = t749 + t752;
  t766 = -1.*t663*t745;
  t767 = -1.*t686*t753;
  t768 = 0.0265*t690;
  t769 = -1.*t560*t686;
  t770 = t769 + t670;
  t775 = -0.0225*t770;
  t776 = t766 + t767 + t768 + t775;
  t748 = -1.*t559*t745;
  t754 = -1.*t663*t753;
  t758 = -1.*t560*t663;
  t759 = t758 + t624;
  t760 = -0.0225*t759;
  t762 = 0.0265*t671;
  t763 = t748 + t754 + t760 + t762;
  t675 = -1.*t506*t492;
  t677 = t443*t542;
  t679 = t675 + t677 + t664 + t665 + t669 + t672;
  t494 = -1.*t443*t492;
  t544 = -1.*t506*t542;
  t580 = -1.*t559*t573;
  t614 = -1.*t593*t607;
  t625 = t623 + t624;
  t626 = -0.0225*t625;
  t631 = t593*t569;
  t641 = t627 + t631;
  t643 = 0.0265*t641;
  t646 = t494 + t544 + t580 + t614 + t626 + t643;
  t796 = 1/Sqrt(t732);
  t824 = Cos(var1[5]);
  t825 = -1.*t824;
  t826 = 1. + t825;
  t828 = Sin(var1[5]);
  t839 = Cos(var1[6]);
  t840 = -1.*t839;
  t841 = 1. + t840;
  t843 = Sin(var1[6]);
  t836 = t506*t824;
  t837 = -1.*t443*t828;
  t838 = t836 + t837;
  t847 = -1.*t824*t443;
  t848 = -1.*t506*t828;
  t849 = t847 + t848;
  t827 = -0.0695*t826;
  t829 = -0.0265*t828;
  t830 = t827 + t829;
  t832 = -0.0265*t826;
  t833 = 0.0695*t828;
  t834 = t832 + t833;
  t842 = -0.2375*t841;
  t844 = -0.0265*t843;
  t845 = t842 + t844;
  t850 = -0.0265*t841;
  t851 = 0.2375*t843;
  t852 = t850 + t851;
  t854 = t839*t838;
  t866 = t824*t443;
  t867 = t506*t828;
  t868 = t866 + t867;
  t831 = -1.*t506*t830;
  t835 = t443*t834;
  t846 = -1.*t838*t845;
  t853 = -1.*t849*t852;
  t855 = t849*t843;
  t856 = t854 + t855;
  t857 = -0.0225*t856;
  t858 = t839*t849;
  t859 = -1.*t838*t843;
  t860 = t858 + t859;
  t861 = 0.0265*t860;
  t862 = t831 + t835 + t846 + t853 + t857 + t861;
  t871 = -1.*t868*t843;
  t872 = t854 + t871;
  t882 = -0.0265*t839;
  t883 = -0.2375*t843;
  t884 = t882 + t883;
  t886 = 0.2375*t839;
  t887 = t886 + t844;
  t864 = -1.*t443*t830;
  t865 = -1.*t506*t834;
  t869 = -1.*t868*t845;
  t870 = -1.*t838*t852;
  t873 = 0.0265*t872;
  t874 = t839*t868;
  t875 = t838*t843;
  t876 = t874 + t875;
  t877 = -0.0225*t876;
  t878 = t864 + t865 + t869 + t870 + t873 + t877;
  t905 = -1.*t506*t824;
  t906 = t443*t828;
  t907 = t905 + t906;
  t898 = -1.*t849*t843;
  t895 = -1.*t838*t884;
  t896 = -1.*t849*t887;
  t897 = -1.*t839*t838;
  t899 = t897 + t898;
  t900 = 0.0265*t899;
  t901 = -0.0225*t860;
  t902 = t895 + t896 + t900 + t901;
  t909 = t839*t907;
  t910 = t909 + t898;
  t885 = -1.*t868*t884;
  t888 = -1.*t838*t887;
  t889 = -0.0225*t872;
  t890 = -1.*t839*t868;
  t891 = t890 + t859;
  t892 = 0.0265*t891;
  t893 = t885 + t888 + t889 + t892;
  t863 = Power(t862,2);
  t879 = Power(t878,2);
  t880 = 0.00085849 + t863 + t879;
  t918 = t443*t830;
  t919 = t506*t834;
  t920 = -1.*t849*t845;
  t921 = -1.*t907*t852;
  t922 = 0.0265*t910;
  t923 = t907*t843;
  t924 = t858 + t923;
  t925 = -0.0225*t924;
  t926 = t918 + t919 + t920 + t921 + t922 + t925;
  t933 = Power(t880,-1.5);
  t938 = -0.0265*t824;
  t939 = -0.0695*t828;
  t940 = t938 + t939;
  t942 = 0.0695*t824;
  t943 = t942 + t829;
  t934 = 2.*t862*t878;
  t935 = 2.*t862*t926;
  t936 = t934 + t935;
  t881 = 1/Sqrt(t880);
  t941 = -1.*t443*t940;
  t944 = -1.*t506*t943;
  t945 = t941 + t944 + t846 + t853 + t857 + t861;
  t913 = -1.*t907*t843;
  t947 = -1.*t506*t940;
  t948 = t443*t943;
  t949 = t947 + t948 + t920 + t921 + t922 + t925;
  t958 = -1.*t907*t845;
  t959 = -1.*t868*t852;
  t960 = t868*t843;
  t961 = t909 + t960;
  t962 = -0.0225*t961;
  t963 = t874 + t913;
  t964 = 0.0265*t963;
  t674 = 2.*t646*t673;
  t697 = 2.*t679*t696;
  t698 = t674 + t697;
  t765 = 2.*t673*t763;
  t777 = 2.*t776*t696;
  t778 = t765 + t777;
  t991 = -0.0265*t458;
  t992 = 0.0695*t482;
  t993 = t991 + t992;
  t794 = 2.*t763*t696;
  t798 = 2.*t724*t679;
  t801 = 2.*t679*t673;
  t804 = 2.*t646*t696;
  t805 = t443*t492;
  t807 = t506*t542;
  t808 = t805 + t807 + t687 + t688 + t691 + t695;
  t810 = 2.*t808*t696;
  t811 = t798 + t801 + t804 + t810;
  t1002 = -1.*t593*t745;
  t1003 = -1.*t559*t753;
  t1004 = 0.0265*t625;
  t1005 = -1.*t560*t559;
  t1006 = t1005 + t720;
  t1007 = -0.0225*t1006;
  t1008 = t1002 + t1003 + t1004 + t1007;
  t1009 = 2.*t1008*t673;
  t1010 = 2.*t646*t763;
  t1011 = 2.*t679*t776;
  t1012 = t1009 + t1010 + t1011 + t794;
  t781 = 2.*t724*t776;
  t782 = 2.*t673*t776;
  t783 = -1.*t686*t745;
  t784 = -1.*t593*t753;
  t785 = -1.*t560*t593;
  t787 = t785 + t689;
  t788 = -0.0225*t787;
  t789 = 0.0265*t721;
  t791 = t783 + t784 + t788 + t789;
  t792 = 2.*t791*t696;
  t795 = t781 + t782 + t792 + t794;
  t1026 = -0.0265*t560;
  t1027 = 0.2375*t569;
  t1028 = t1026 + t1027;
  t903 = 2.*t902*t878;
  t904 = -1.*t849*t884;
  t908 = -1.*t907*t887;
  t911 = -0.0225*t910;
  t912 = -1.*t839*t849;
  t914 = t912 + t913;
  t915 = 0.0265*t914;
  t916 = t904 + t908 + t911 + t915;
  t917 = 2.*t862*t916;
  t930 = 2.*t902*t862;
  t931 = 2.*t893*t878;
  t932 = t930 + t931;
  t946 = 2.*t945*t878;
  t950 = 2.*t862*t949;
  t951 = t946 + t950;
  t956 = t443*t940;
  t1053 = -0.0695*t824;
  t1054 = 0.0265*t828;
  t1055 = t1053 + t1054;
  t955 = 2.*t945*t862;
  t957 = t506*t943;
  t965 = t956 + t957 + t958 + t959 + t962 + t964;
  t966 = 2.*t862*t965;
  t967 = 2.*t878*t949;
  t968 = 2.*t949*t926;
  t969 = t955 + t966 + t967 + t968;
  t1075 = -0.2375*t839;
  t1076 = 0.0265*t843;
  t1077 = t1075 + t1076;
  t1044 = 2.*t945*t893;
  t1045 = 2.*t902*t949;
  t1046 = t1044 + t903 + t917 + t1045;
  t894 = 2.*t862*t893;
  t927 = 2.*t902*t926;
  t928 = t894 + t903 + t917 + t927;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=-0.25*Power(t728,2)*t735*var1[9] + 0.5*(2.*(t580 + t614 + t626 + t643 - 1.*t443*t652 + t506*t658)*t696 + 2.*t673*t724 + 2.*Power(t724,2) + 2.*t731)*t796*var1[9] - 0.25*t698*t728*t735*var1[10] + 0.5*t796*t811*var1[10] - 0.25*t728*t735*t778*var1[11] + 0.5*t795*t796*var1[11];
  p_output1[5]=-0.25*t933*Power(t936,2)*var1[9] + 0.5*t881*(2.*t863 + 2.*t878*t926 + 2.*Power(t926,2) + 2.*t862*(t506*t830 - 1.*t443*t834 + t958 + t959 + t962 + t964))*var1[9] - 0.25*t933*t936*t951*var1[12] + 0.5*t881*t969*var1[12] + 0.5*t881*t928*var1[13] - 0.25*t932*t933*t936*var1[13];
  p_output1[6]=-0.25*t698*t728*t735*var1[9] + 0.5*t796*t811*var1[9] - 0.25*Power(t698,2)*t735*var1[10] + 0.5*t796*(2.*Power(t646,2) + 2.*Power(t679,2) + 2.*t696*(t544 + t580 + t614 + t626 + t643 + t443*t993) + 2.*t673*(-1.*t443*t542 + t702 + t705 + t714 + t723 - 1.*t506*t993))*var1[10] - 0.25*t698*t735*t778*var1[11] + 0.5*t1012*t796*var1[11];
  p_output1[7]=0;
  p_output1[8]=-0.25*t728*t735*t778*var1[9] + 0.5*t795*t796*var1[9] - 0.25*t698*t735*t778*var1[10] + 0.5*t1012*t796*var1[10] - 0.25*t735*Power(t778,2)*var1[11] + 0.5*(2.*t673*(t1003 - 1.*t1028*t663 - 0.0225*(t1005 + t693) + 0.0265*t759) + 2.*Power(t763,2) + 2.*t696*(-1.*t1028*t686 + t754 - 0.0225*(t710 + t758) + 0.0265*t770) + 2.*Power(t776,2))*t796*var1[11];
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=-0.25*t933*t936*t951*var1[9] + 0.5*t881*t969*var1[9] - 0.25*t933*Power(t951,2)*var1[12] + 0.5*t881*(2.*Power(t945,2) + 2.*t878*(-1.*t1055*t443 + t920 + t921 + t922 + t925 + t947) + 2.*Power(t949,2) + 2.*t862*(-1.*t1055*t506 + t956 + t958 + t959 + t962 + t964))*var1[12] + 0.5*t1046*t881*var1[13] - 0.25*t932*t933*t951*var1[13];
  p_output1[12]=0;
  p_output1[13]=0.5*t881*t928*var1[9] - 0.25*t932*t933*t936*var1[9] + 0.5*t1046*t881*var1[12] - 0.25*t932*t933*t951*var1[12] - 0.25*Power(t932,2)*t933*var1[13] + 0.5*t881*(2.*Power(t893,2) + 2.*Power(t902,2) + 2.*t862*(-1.*t1077*t838 - 0.0225*t899 + t904 + 0.0265*(t875 + t912)) + 2.*t878*(-1.*t1077*t868 - 0.0225*t891 + t895 + 0.0265*(t897 + t960)))*var1[13];
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

#include "dJ_LegL.hh"

namespace SymFunction
{

void dJ_LegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
