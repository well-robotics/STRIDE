/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:53 GMT-05:00
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
  double t2356;
  double t2392;
  double t780;
  double t2441;
  double t2473;
  double t2474;
  double t2476;
  double t2483;
  double t2505;
  double t2507;
  double t2517;
  double t2463;
  double t2466;
  double t2470;
  double t2412;
  double t2559;
  double t2560;
  double t2479;
  double t2487;
  double t2495;
  double t2528;
  double t2530;
  double t2532;
  double t2575;
  double t2576;
  double t2582;
  double t2549;
  double t2374;
  double t2413;
  double t2444;
  double t2446;
  double t2451;
  double t2584;
  double t2590;
  double t2591;
  double t2594;
  double t2595;
  double t2596;
  double t2598;
  double t2600;
  double t2603;
  double t2561;
  double t2562;
  double t2567;
  double t2569;
  double t2572;
  double t2620;
  double t2622;
  double t2623;
  double t2563;
  double t2573;
  double t2606;
  double t2616;
  double t2617;
  double t2626;
  double t2627;
  double t2629;
  double t2634;
  double t2637;
  double t2641;
  double t2643;
  double t2644;
  double t2646;
  double t2647;
  double t2537;
  double t2651;
  double t2653;
  double t2654;
  double t2655;
  double t2539;
  double t2685;
  double t2687;
  double t2689;
  double t2691;
  double t2692;
  double t2716;
  double t2720;
  double t2723;
  double t2725;
  double t2719;
  double t2721;
  double t2722;
  double t2729;
  double t2731;
  double t2732;
  double t2747;
  double t2748;
  double t2735;
  double t2758;
  double t2761;
  double t2742;
  double t2743;
  double t2744;
  double t2724;
  double t2726;
  double t2727;
  double t2734;
  double t2736;
  double t2771;
  double t2782;
  double t2783;
  double t2784;
  double t2749;
  double t2750;
  double t2751;
  double t2753;
  double t2754;
  double t2756;
  double t2762;
  double t2763;
  double t2766;
  double t2767;
  double t2768;
  double t2787;
  double t2788;
  double t2752;
  double t2757;
  double t2764;
  double t2770;
  double t2774;
  double t2775;
  double t2776;
  double t2777;
  double t2778;
  double t2795;
  double t2796;
  double t2797;
  double t2798;
  double t2799;
  double t2800;
  double t2801;
  double t2804;
  double t2805;
  double t2806;
  double t2809;
  double t2810;
  double t2811;
  double t2812;
  double t2814;
  double t2815;
  double t2816;
  double t2818;
  double t2819;
  double t2826;
  double t2827;
  double t2828;
  double t2739;
  double t2825;
  double t2829;
  double t2830;
  double t2831;
  double t2832;
  double t2833;
  double t2834;
  double t2835;
  t2356 = Cos(var1[3]);
  t2392 = Sin(var1[3]);
  t780 = Sin(var1[2]);
  t2441 = Cos(var1[2]);
  t2473 = Cos(var1[4]);
  t2474 = -1.*t2473;
  t2476 = 1. + t2474;
  t2483 = Sin(var1[4]);
  t2505 = -1.*t2441*t2356;
  t2507 = -1.*t780*t2392;
  t2517 = t2505 + t2507;
  t2463 = t2356*t780;
  t2466 = -1.*t2441*t2392;
  t2470 = t2463 + t2466;
  t2412 = -0.0695*t2392;
  t2559 = -1.*t2356;
  t2560 = 1. + t2559;
  t2479 = -0.0265*t2476;
  t2487 = -0.2375*t2483;
  t2495 = t2479 + t2487;
  t2528 = -0.2375*t2476;
  t2530 = 0.0265*t2483;
  t2532 = t2528 + t2530;
  t2575 = t2441*t2356;
  t2576 = t780*t2392;
  t2582 = t2575 + t2576;
  t2549 = t2473*t2470;
  t2374 = 0.0265*t2356;
  t2413 = t2374 + t2412;
  t2444 = -0.0695*t2356;
  t2446 = -0.0265*t2392;
  t2451 = t2444 + t2446;
  t2584 = -1.*t2582*t2495;
  t2590 = -1.*t2470*t2532;
  t2591 = t2473*t2582;
  t2594 = t2470*t2483;
  t2595 = t2591 + t2594;
  t2596 = 0.0265*t2595;
  t2598 = -1.*t2582*t2483;
  t2600 = t2549 + t2598;
  t2603 = -0.0325*t2600;
  t2561 = -0.0265*t2560;
  t2562 = t2561 + t2412;
  t2567 = -0.0695*t2560;
  t2569 = 0.0265*t2392;
  t2572 = t2567 + t2569;
  t2620 = -1.*t2356*t780;
  t2622 = t2441*t2392;
  t2623 = t2620 + t2622;
  t2563 = -1.*t2441*t2562;
  t2573 = -1.*t780*t2572;
  t2606 = t2563 + t2573 + t2584 + t2590 + t2596 + t2603;
  t2616 = t780*t2562;
  t2617 = -1.*t2441*t2572;
  t2626 = -1.*t2623*t2495;
  t2627 = -1.*t2582*t2532;
  t2629 = -1.*t2623*t2483;
  t2634 = t2591 + t2629;
  t2637 = -0.0325*t2634;
  t2641 = t2473*t2623;
  t2643 = t2582*t2483;
  t2644 = t2641 + t2643;
  t2646 = 0.0265*t2644;
  t2647 = t2616 + t2617 + t2626 + t2627 + t2637 + t2646;
  t2537 = t2473*t2517;
  t2651 = Power(t2606,2);
  t2653 = Power(t2647,2);
  t2654 = 0.00085849 + t2651 + t2653;
  t2655 = 1/Sqrt(t2654);
  t2539 = -1.*t2470*t2483;
  t2685 = 0.0265*t2473;
  t2687 = t2685 + t2487;
  t2689 = -0.2375*t2473;
  t2691 = -0.0265*t2483;
  t2692 = t2689 + t2691;
  t2716 = Cos(var1[5]);
  t2720 = Sin(var1[5]);
  t2723 = Cos(var1[6]);
  t2725 = Sin(var1[6]);
  t2719 = t2441*t2716;
  t2721 = -1.*t780*t2720;
  t2722 = t2719 + t2721;
  t2729 = -1.*t2716*t780;
  t2731 = -1.*t2441*t2720;
  t2732 = t2729 + t2731;
  t2747 = -1.*t2716;
  t2748 = 1. + t2747;
  t2735 = -0.0265*t2725;
  t2758 = -1.*t2723;
  t2761 = 1. + t2758;
  t2742 = t2723*t2732;
  t2743 = -1.*t2722*t2725;
  t2744 = t2742 + t2743;
  t2724 = -0.0265*t2723;
  t2726 = -0.2375*t2725;
  t2727 = t2724 + t2726;
  t2734 = 0.2375*t2723;
  t2736 = t2734 + t2735;
  t2771 = t2723*t2722;
  t2782 = t2716*t780;
  t2783 = t2441*t2720;
  t2784 = t2782 + t2783;
  t2749 = -0.0695*t2748;
  t2750 = -0.0265*t2720;
  t2751 = t2749 + t2750;
  t2753 = -0.0265*t2748;
  t2754 = 0.0695*t2720;
  t2756 = t2753 + t2754;
  t2762 = -0.2375*t2761;
  t2763 = t2762 + t2735;
  t2766 = -0.0265*t2761;
  t2767 = 0.2375*t2725;
  t2768 = t2766 + t2767;
  t2787 = -1.*t2784*t2725;
  t2788 = t2771 + t2787;
  t2752 = -1.*t2441*t2751;
  t2757 = t780*t2756;
  t2764 = -1.*t2722*t2763;
  t2770 = -1.*t2732*t2768;
  t2774 = t2732*t2725;
  t2775 = t2771 + t2774;
  t2776 = -0.0325*t2775;
  t2777 = 0.0265*t2744;
  t2778 = t2752 + t2757 + t2764 + t2770 + t2776 + t2777;
  t2795 = -1.*t780*t2751;
  t2796 = -1.*t2441*t2756;
  t2797 = -1.*t2784*t2763;
  t2798 = -1.*t2722*t2768;
  t2799 = 0.0265*t2788;
  t2800 = t2723*t2784;
  t2801 = t2722*t2725;
  t2804 = t2800 + t2801;
  t2805 = -0.0325*t2804;
  t2806 = t2795 + t2796 + t2797 + t2798 + t2799 + t2805;
  t2809 = Power(t2778,2);
  t2810 = Power(t2806,2);
  t2811 = 0.00085849 + t2809 + t2810;
  t2812 = 1/Sqrt(t2811);
  t2814 = -0.0265*t2716;
  t2815 = -0.0695*t2720;
  t2816 = t2814 + t2815;
  t2818 = 0.0695*t2716;
  t2819 = t2818 + t2750;
  t2826 = -1.*t2441*t2716;
  t2827 = t780*t2720;
  t2828 = t2826 + t2827;
  t2739 = -1.*t2732*t2725;
  t2825 = -1.*t2732*t2763;
  t2829 = -1.*t2828*t2768;
  t2830 = t2723*t2828;
  t2831 = t2830 + t2739;
  t2832 = 0.0265*t2831;
  t2833 = t2828*t2725;
  t2834 = t2742 + t2833;
  t2835 = -0.0325*t2834;
  p_output1[0]=0.5*t2655*(2.*t2606*t2647 + 2.*t2647*(-1.*t2495*t2517 + t2441*t2562 - 1.*t2532*t2623 + 0.0265*(t2537 + t2483*t2623) - 0.0325*(-1.*t2483*t2517 + t2641) + t2572*t780))*var1[9] + 0.5*t2655*(2.*t2606*(-1.*t2441*t2451 - 1.*t2470*t2495 - 1.*t2517*t2532 - 0.0325*(t2537 + t2539) + 0.0265*(t2483*t2517 + t2549) - 1.*t2413*t780) + 2.*t2647*(-1.*t2413*t2441 + t2584 + t2590 + t2596 + t2603 + t2451*t780))*var1[10] + 0.5*t2655*(2.*t2606*(-0.0325*(t2539 - 1.*t2473*t2582) + 0.0265*t2600 - 1.*t2470*t2687 - 1.*t2582*t2692) + 2.*t2647*(-0.0325*(t2598 - 1.*t2473*t2623) + 0.0265*t2634 - 1.*t2582*t2687 - 1.*t2623*t2692))*var1[11];
  p_output1[1]=0.5*t2812*(2.*t2778*t2806 + 2.*t2778*(t2441*t2756 + t2825 + t2829 + t2832 + t2835 + t2751*t780))*var1[9] + 0.5*t2812*(2.*t2806*(t2764 + t2770 + t2776 + t2777 - 1.*t2441*t2819 - 1.*t2816*t780) + 2.*t2778*(-1.*t2441*t2816 + t2825 + t2829 + t2832 + t2835 + t2819*t780))*var1[12] + 0.5*(2.*(-1.*t2722*t2727 - 1.*t2732*t2736 + 0.0265*(-1.*t2722*t2723 + t2739) - 0.0325*t2744)*t2778 + 2.*(-1.*t2722*t2736 - 1.*t2727*t2784 + 0.0265*(t2743 - 1.*t2723*t2784) - 0.0325*t2788)*t2806)*t2812*var1[13];
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
