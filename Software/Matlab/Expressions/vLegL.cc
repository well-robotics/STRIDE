/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:06 GMT-05:00
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
static void output1(double *p_output1,const double *var1)
{
  double t2428;
  double t2464;
  double t2409;
  double t2499;
  double t2530;
  double t2532;
  double t2534;
  double t2543;
  double t2563;
  double t2573;
  double t2584;
  double t2522;
  double t2526;
  double t2527;
  double t2469;
  double t2616;
  double t2617;
  double t2539;
  double t2551;
  double t2553;
  double t2585;
  double t2588;
  double t2591;
  double t2632;
  double t2638;
  double t2640;
  double t2608;
  double t2445;
  double t2483;
  double t2502;
  double t2507;
  double t2515;
  double t2646;
  double t2647;
  double t2650;
  double t2651;
  double t2652;
  double t2654;
  double t2656;
  double t2659;
  double t2662;
  double t2618;
  double t2619;
  double t2625;
  double t2628;
  double t2629;
  double t2678;
  double t2679;
  double t2682;
  double t2623;
  double t2631;
  double t2663;
  double t2673;
  double t2674;
  double t2683;
  double t2684;
  double t2690;
  double t2693;
  double t2695;
  double t2699;
  double t2700;
  double t2702;
  double t2703;
  double t2704;
  double t2595;
  double t2709;
  double t2710;
  double t2711;
  double t2713;
  double t2600;
  double t2743;
  double t2744;
  double t2747;
  double t2748;
  double t2749;
  double t2775;
  double t2777;
  double t2780;
  double t2782;
  double t2776;
  double t2778;
  double t2779;
  double t2787;
  double t2788;
  double t2790;
  double t2804;
  double t2805;
  double t2792;
  double t2815;
  double t2818;
  double t2799;
  double t2800;
  double t2801;
  double t2781;
  double t2783;
  double t2784;
  double t2791;
  double t2793;
  double t2828;
  double t2839;
  double t2840;
  double t2841;
  double t2806;
  double t2807;
  double t2808;
  double t2810;
  double t2811;
  double t2813;
  double t2819;
  double t2820;
  double t2823;
  double t2824;
  double t2825;
  double t2844;
  double t2845;
  double t2809;
  double t2814;
  double t2821;
  double t2827;
  double t2831;
  double t2832;
  double t2833;
  double t2834;
  double t2835;
  double t2852;
  double t2853;
  double t2854;
  double t2855;
  double t2856;
  double t2857;
  double t2858;
  double t2861;
  double t2862;
  double t2863;
  double t2866;
  double t2867;
  double t2868;
  double t2869;
  double t2871;
  double t2872;
  double t2873;
  double t2875;
  double t2876;
  double t2883;
  double t2884;
  double t2885;
  double t2796;
  double t2882;
  double t2886;
  double t2887;
  double t2888;
  double t2889;
  double t2890;
  double t2891;
  double t2892;
  t2428 = Cos(var1[3]);
  t2464 = Sin(var1[3]);
  t2409 = Sin(var1[2]);
  t2499 = Cos(var1[2]);
  t2530 = Cos(var1[4]);
  t2532 = -1.*t2530;
  t2534 = 1. + t2532;
  t2543 = Sin(var1[4]);
  t2563 = -1.*t2499*t2428;
  t2573 = -1.*t2409*t2464;
  t2584 = t2563 + t2573;
  t2522 = t2428*t2409;
  t2526 = -1.*t2499*t2464;
  t2527 = t2522 + t2526;
  t2469 = -0.0695*t2464;
  t2616 = -1.*t2428;
  t2617 = 1. + t2616;
  t2539 = -0.0265*t2534;
  t2551 = -0.2375*t2543;
  t2553 = t2539 + t2551;
  t2585 = -0.2375*t2534;
  t2588 = 0.0265*t2543;
  t2591 = t2585 + t2588;
  t2632 = t2499*t2428;
  t2638 = t2409*t2464;
  t2640 = t2632 + t2638;
  t2608 = t2530*t2527;
  t2445 = 0.0265*t2428;
  t2483 = t2445 + t2469;
  t2502 = -0.0695*t2428;
  t2507 = -0.0265*t2464;
  t2515 = t2502 + t2507;
  t2646 = -1.*t2640*t2553;
  t2647 = -1.*t2527*t2591;
  t2650 = t2530*t2640;
  t2651 = t2527*t2543;
  t2652 = t2650 + t2651;
  t2654 = 0.0265*t2652;
  t2656 = -1.*t2640*t2543;
  t2659 = t2608 + t2656;
  t2662 = 0.0115*t2659;
  t2618 = -0.0265*t2617;
  t2619 = t2618 + t2469;
  t2625 = -0.0695*t2617;
  t2628 = 0.0265*t2464;
  t2629 = t2625 + t2628;
  t2678 = -1.*t2428*t2409;
  t2679 = t2499*t2464;
  t2682 = t2678 + t2679;
  t2623 = -1.*t2499*t2619;
  t2631 = -1.*t2409*t2629;
  t2663 = t2623 + t2631 + t2646 + t2647 + t2654 + t2662;
  t2673 = t2409*t2619;
  t2674 = -1.*t2499*t2629;
  t2683 = -1.*t2682*t2553;
  t2684 = -1.*t2640*t2591;
  t2690 = -1.*t2682*t2543;
  t2693 = t2650 + t2690;
  t2695 = 0.0115*t2693;
  t2699 = t2530*t2682;
  t2700 = t2640*t2543;
  t2702 = t2699 + t2700;
  t2703 = 0.0265*t2702;
  t2704 = t2673 + t2674 + t2683 + t2684 + t2695 + t2703;
  t2595 = t2530*t2584;
  t2709 = Power(t2663,2);
  t2710 = Power(t2704,2);
  t2711 = 0.00085849 + t2709 + t2710;
  t2713 = 1/Sqrt(t2711);
  t2600 = -1.*t2527*t2543;
  t2743 = 0.0265*t2530;
  t2744 = t2743 + t2551;
  t2747 = -0.2375*t2530;
  t2748 = -0.0265*t2543;
  t2749 = t2747 + t2748;
  t2775 = Cos(var1[5]);
  t2777 = Sin(var1[5]);
  t2780 = Cos(var1[6]);
  t2782 = Sin(var1[6]);
  t2776 = t2499*t2775;
  t2778 = -1.*t2409*t2777;
  t2779 = t2776 + t2778;
  t2787 = -1.*t2775*t2409;
  t2788 = -1.*t2499*t2777;
  t2790 = t2787 + t2788;
  t2804 = -1.*t2775;
  t2805 = 1. + t2804;
  t2792 = -0.0265*t2782;
  t2815 = -1.*t2780;
  t2818 = 1. + t2815;
  t2799 = t2780*t2790;
  t2800 = -1.*t2779*t2782;
  t2801 = t2799 + t2800;
  t2781 = -0.0265*t2780;
  t2783 = -0.2375*t2782;
  t2784 = t2781 + t2783;
  t2791 = 0.2375*t2780;
  t2793 = t2791 + t2792;
  t2828 = t2780*t2779;
  t2839 = t2775*t2409;
  t2840 = t2499*t2777;
  t2841 = t2839 + t2840;
  t2806 = -0.0695*t2805;
  t2807 = -0.0265*t2777;
  t2808 = t2806 + t2807;
  t2810 = -0.0265*t2805;
  t2811 = 0.0695*t2777;
  t2813 = t2810 + t2811;
  t2819 = -0.2375*t2818;
  t2820 = t2819 + t2792;
  t2823 = -0.0265*t2818;
  t2824 = 0.2375*t2782;
  t2825 = t2823 + t2824;
  t2844 = -1.*t2841*t2782;
  t2845 = t2828 + t2844;
  t2809 = -1.*t2499*t2808;
  t2814 = t2409*t2813;
  t2821 = -1.*t2779*t2820;
  t2827 = -1.*t2790*t2825;
  t2831 = t2790*t2782;
  t2832 = t2828 + t2831;
  t2833 = 0.0115*t2832;
  t2834 = 0.0265*t2801;
  t2835 = t2809 + t2814 + t2821 + t2827 + t2833 + t2834;
  t2852 = -1.*t2409*t2808;
  t2853 = -1.*t2499*t2813;
  t2854 = -1.*t2841*t2820;
  t2855 = -1.*t2779*t2825;
  t2856 = 0.0265*t2845;
  t2857 = t2780*t2841;
  t2858 = t2779*t2782;
  t2861 = t2857 + t2858;
  t2862 = 0.0115*t2861;
  t2863 = t2852 + t2853 + t2854 + t2855 + t2856 + t2862;
  t2866 = Power(t2835,2);
  t2867 = Power(t2863,2);
  t2868 = 0.00085849 + t2866 + t2867;
  t2869 = 1/Sqrt(t2868);
  t2871 = -0.0265*t2775;
  t2872 = -0.0695*t2777;
  t2873 = t2871 + t2872;
  t2875 = 0.0695*t2775;
  t2876 = t2875 + t2807;
  t2883 = -1.*t2499*t2775;
  t2884 = t2409*t2777;
  t2885 = t2883 + t2884;
  t2796 = -1.*t2790*t2782;
  t2882 = -1.*t2790*t2820;
  t2886 = -1.*t2885*t2825;
  t2887 = t2780*t2885;
  t2888 = t2887 + t2796;
  t2889 = 0.0265*t2888;
  t2890 = t2885*t2782;
  t2891 = t2799 + t2890;
  t2892 = 0.0115*t2891;
  p_output1[0]=0.5*(2.*t2663*t2704 + 2.*(-1.*t2553*t2584 + t2499*t2619 + t2409*t2629 - 1.*t2591*t2682 + 0.0265*(t2595 + t2543*t2682) + 0.0115*(-1.*t2543*t2584 + t2699))*t2704)*t2713*var1[9] + 0.5*(2.*(-1.*t2409*t2483 - 1.*t2499*t2515 - 1.*t2527*t2553 - 1.*t2584*t2591 + 0.0115*(t2595 + t2600) + 0.0265*(t2543*t2584 + t2608))*t2663 + 2.*(-1.*t2483*t2499 + t2409*t2515 + t2646 + t2647 + t2654 + t2662)*t2704)*t2713*var1[10] + 0.5*t2713*(2.*t2663*(0.0115*(t2600 - 1.*t2530*t2640) + 0.0265*t2659 - 1.*t2527*t2744 - 1.*t2640*t2749) + 2.*t2704*(0.0115*(t2656 - 1.*t2530*t2682) + 0.0265*t2693 - 1.*t2640*t2744 - 1.*t2682*t2749))*var1[11];
  p_output1[1]=0.5*t2869*(2.*t2835*t2863 + 2.*t2835*(t2409*t2808 + t2499*t2813 + t2882 + t2886 + t2889 + t2892))*var1[9] + 0.5*t2869*(2.*t2863*(t2821 + t2827 + t2833 + t2834 - 1.*t2409*t2873 - 1.*t2499*t2876) + 2.*t2835*(-1.*t2499*t2873 + t2409*t2876 + t2882 + t2886 + t2889 + t2892))*var1[12] + 0.5*(2.*(-1.*t2779*t2784 - 1.*t2790*t2793 + 0.0265*(-1.*t2779*t2780 + t2796) + 0.0115*t2801)*t2835 + 2.*(-1.*t2779*t2793 - 1.*t2784*t2841 + 0.0265*(t2800 - 1.*t2780*t2841) + 0.0115*t2845)*t2863)*t2869*var1[13];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "vLegL.hh"

namespace SymFunction
{

void vLegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
