/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:03 GMT-05:00
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
  double t2365;
  double t2390;
  double t2406;
  double t2424;
  double t2342;
  double t2451;
  double t2492;
  double t2494;
  double t2495;
  double t2498;
  double t2485;
  double t2488;
  double t2490;
  double t2503;
  double t2507;
  double t2512;
  double t2409;
  double t2428;
  double t2431;
  double t2464;
  double t2469;
  double t2480;
  double t2497;
  double t2499;
  double t2501;
  double t2515;
  double t2519;
  double t2520;
  double t2561;
  double t2563;
  double t2564;
  double t2533;
  double t2573;
  double t2600;
  double t2601;
  double t2602;
  double t2553;
  double t2556;
  double t2557;
  double t2571;
  double t2580;
  double t2584;
  double t2585;
  double t2588;
  double t2591;
  double t2593;
  double t2594;
  double t2597;
  double t2598;
  double t2599;
  double t2603;
  double t2604;
  double t2605;
  double t2607;
  double t2608;
  double t2609;
  double t2610;
  double t2611;
  double t2612;
  double t2620;
  double t2621;
  double t2622;
  double t2624;
  double t2635;
  double t2636;
  double t2637;
  double t2639;
  double t2632;
  double t2633;
  double t2634;
  double t2643;
  double t2644;
  double t2645;
  double t2623;
  double t2625;
  double t2626;
  double t2628;
  double t2629;
  double t2630;
  double t2638;
  double t2640;
  double t2641;
  double t2646;
  double t2647;
  double t2648;
  double t2650;
  double t2662;
  double t2663;
  double t2664;
  double t2627;
  double t2631;
  double t2642;
  double t2649;
  double t2651;
  double t2652;
  double t2653;
  double t2654;
  double t2655;
  double t2656;
  double t2657;
  double t2658;
  double t2660;
  double t2661;
  double t2665;
  double t2666;
  double t2667;
  double t2668;
  double t2669;
  double t2670;
  double t2672;
  double t2673;
  double t2674;
  double t2677;
  double t2687;
  double t2690;
  double t2691;
  double t2526;
  double t2706;
  double t2708;
  double t2710;
  double t2711;
  double t2712;
  double t2615;
  double t2616;
  double t2617;
  double t2618;
  double t2720;
  double t2741;
  double t2742;
  double t2744;
  double t2745;
  double t2746;
  double t2659;
  double t2678;
  double t2679;
  double t2680;
  double t2767;
  double t2768;
  double t2769;
  double t2772;
  double t2773;
  double t2686;
  double t2692;
  double t2693;
  double t2694;
  double t2695;
  double t2698;
  double t2699;
  double t2700;
  double t2701;
  double t2784;
  double t2785;
  double t2786;
  double t2788;
  double t2789;
  t2365 = Cos(var1[3]);
  t2390 = -1.*t2365;
  t2406 = 1. + t2390;
  t2424 = Sin(var1[3]);
  t2342 = Cos(var1[2]);
  t2451 = Sin(var1[2]);
  t2492 = Cos(var1[4]);
  t2494 = -1.*t2492;
  t2495 = 1. + t2494;
  t2498 = Sin(var1[4]);
  t2485 = -1.*t2342*t2365;
  t2488 = -1.*t2451*t2424;
  t2490 = t2485 + t2488;
  t2503 = -1.*t2365*t2451;
  t2507 = t2342*t2424;
  t2512 = t2503 + t2507;
  t2409 = -0.0265*t2406;
  t2428 = -0.0695*t2424;
  t2431 = t2409 + t2428;
  t2464 = -0.0695*t2406;
  t2469 = 0.0265*t2424;
  t2480 = t2464 + t2469;
  t2497 = -0.0265*t2495;
  t2499 = -0.2375*t2498;
  t2501 = t2497 + t2499;
  t2515 = -0.2375*t2495;
  t2519 = 0.0265*t2498;
  t2520 = t2515 + t2519;
  t2561 = t2342*t2365;
  t2563 = t2451*t2424;
  t2564 = t2561 + t2563;
  t2533 = t2492*t2512;
  t2573 = t2492*t2564;
  t2600 = t2365*t2451;
  t2601 = -1.*t2342*t2424;
  t2602 = t2600 + t2601;
  t2553 = t2451*t2431;
  t2556 = -1.*t2342*t2480;
  t2557 = -1.*t2512*t2501;
  t2571 = -1.*t2564*t2520;
  t2580 = -1.*t2512*t2498;
  t2584 = t2573 + t2580;
  t2585 = 0.0115*t2584;
  t2588 = t2564*t2498;
  t2591 = t2533 + t2588;
  t2593 = 0.0265*t2591;
  t2594 = t2553 + t2556 + t2557 + t2571 + t2585 + t2593;
  t2597 = -1.*t2342*t2431;
  t2598 = -1.*t2451*t2480;
  t2599 = -1.*t2564*t2501;
  t2603 = -1.*t2602*t2520;
  t2604 = t2602*t2498;
  t2605 = t2573 + t2604;
  t2607 = 0.0265*t2605;
  t2608 = t2492*t2602;
  t2609 = -1.*t2564*t2498;
  t2610 = t2608 + t2609;
  t2611 = 0.0115*t2610;
  t2612 = t2597 + t2598 + t2599 + t2603 + t2607 + t2611;
  t2620 = Cos(var1[5]);
  t2621 = -1.*t2620;
  t2622 = 1. + t2621;
  t2624 = Sin(var1[5]);
  t2635 = Cos(var1[6]);
  t2636 = -1.*t2635;
  t2637 = 1. + t2636;
  t2639 = Sin(var1[6]);
  t2632 = t2342*t2620;
  t2633 = -1.*t2451*t2624;
  t2634 = t2632 + t2633;
  t2643 = -1.*t2620*t2451;
  t2644 = -1.*t2342*t2624;
  t2645 = t2643 + t2644;
  t2623 = -0.0695*t2622;
  t2625 = -0.0265*t2624;
  t2626 = t2623 + t2625;
  t2628 = -0.0265*t2622;
  t2629 = 0.0695*t2624;
  t2630 = t2628 + t2629;
  t2638 = -0.2375*t2637;
  t2640 = -0.0265*t2639;
  t2641 = t2638 + t2640;
  t2646 = -0.0265*t2637;
  t2647 = 0.2375*t2639;
  t2648 = t2646 + t2647;
  t2650 = t2635*t2634;
  t2662 = t2620*t2451;
  t2663 = t2342*t2624;
  t2664 = t2662 + t2663;
  t2627 = -1.*t2342*t2626;
  t2631 = t2451*t2630;
  t2642 = -1.*t2634*t2641;
  t2649 = -1.*t2645*t2648;
  t2651 = t2645*t2639;
  t2652 = t2650 + t2651;
  t2653 = 0.0115*t2652;
  t2654 = t2635*t2645;
  t2655 = -1.*t2634*t2639;
  t2656 = t2654 + t2655;
  t2657 = 0.0265*t2656;
  t2658 = t2627 + t2631 + t2642 + t2649 + t2653 + t2657;
  t2660 = -1.*t2451*t2626;
  t2661 = -1.*t2342*t2630;
  t2665 = -1.*t2664*t2641;
  t2666 = -1.*t2634*t2648;
  t2667 = -1.*t2664*t2639;
  t2668 = t2650 + t2667;
  t2669 = 0.0265*t2668;
  t2670 = t2635*t2664;
  t2672 = t2634*t2639;
  t2673 = t2670 + t2672;
  t2674 = 0.0115*t2673;
  t2677 = t2660 + t2661 + t2665 + t2666 + t2669 + t2674;
  t2687 = -1.*t2342*t2620;
  t2690 = t2451*t2624;
  t2691 = t2687 + t2690;
  t2526 = t2492*t2490;
  t2706 = 0.0265*t2365;
  t2708 = t2706 + t2428;
  t2710 = -0.0695*t2365;
  t2711 = -0.0265*t2424;
  t2712 = t2710 + t2711;
  t2615 = Power(t2612,2);
  t2616 = Power(t2594,2);
  t2617 = 0.00085849 + t2615 + t2616;
  t2618 = 1/Sqrt(t2617);
  t2720 = -1.*t2602*t2498;
  t2741 = 0.0265*t2492;
  t2742 = t2741 + t2499;
  t2744 = -0.2375*t2492;
  t2745 = -0.0265*t2498;
  t2746 = t2744 + t2745;
  t2659 = Power(t2658,2);
  t2678 = Power(t2677,2);
  t2679 = 0.00085849 + t2659 + t2678;
  t2680 = 1/Sqrt(t2679);
  t2767 = -0.0265*t2620;
  t2768 = -0.0695*t2624;
  t2769 = t2767 + t2768;
  t2772 = 0.0695*t2620;
  t2773 = t2772 + t2625;
  t2686 = -1.*t2645*t2641;
  t2692 = -1.*t2691*t2648;
  t2693 = t2635*t2691;
  t2694 = -1.*t2645*t2639;
  t2695 = t2693 + t2694;
  t2698 = 0.0265*t2695;
  t2699 = t2691*t2639;
  t2700 = t2654 + t2699;
  t2701 = 0.0115*t2700;
  t2784 = -0.0265*t2635;
  t2785 = -0.2375*t2639;
  t2786 = t2784 + t2785;
  t2788 = 0.2375*t2635;
  t2789 = t2788 + t2640;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0.5*(2.*(t2342*t2431 + t2451*t2480 - 1.*t2490*t2501 - 1.*t2512*t2520 + 0.0265*(t2498*t2512 + t2526) + 0.0115*(-1.*t2490*t2498 + t2533))*t2594 + 2.*t2594*t2612)*t2618;
  p_output1[5]=0.5*t2680*(2.*t2658*t2677 + 2.*t2658*(t2451*t2626 + t2342*t2630 + t2686 + t2692 + t2698 + t2701));
  p_output1[6]=0.5*t2618*(2.*t2594*(t2599 + t2603 + t2607 + t2611 - 1.*t2342*t2708 + t2451*t2712) + 2.*t2612*(-1.*t2490*t2520 - 1.*t2501*t2602 + 0.0265*(t2490*t2498 + t2608) - 1.*t2451*t2708 - 1.*t2342*t2712 + 0.0115*(t2526 + t2720)));
  p_output1[7]=0;
  p_output1[8]=0.5*t2618*(2.*t2594*(0.0265*t2584 + 0.0115*(-1.*t2492*t2512 + t2609) - 1.*t2564*t2742 - 1.*t2512*t2746) + 2.*t2612*(0.0265*t2610 + 0.0115*(-1.*t2492*t2564 + t2720) - 1.*t2602*t2742 - 1.*t2564*t2746));
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=0.5*t2680*(2.*t2677*(t2642 + t2649 + t2653 + t2657 - 1.*t2451*t2769 - 1.*t2342*t2773) + 2.*t2658*(t2686 + t2692 + t2698 + t2701 - 1.*t2342*t2769 + t2451*t2773));
  p_output1[12]=0;
  p_output1[13]=0.5*t2680*(2.*t2677*(0.0265*(t2655 - 1.*t2635*t2664) + 0.0115*t2668 - 1.*t2664*t2786 - 1.*t2634*t2789) + 2.*t2658*(0.0115*t2656 + 0.0265*(-1.*t2634*t2635 + t2694) - 1.*t2634*t2786 - 1.*t2645*t2789));
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
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
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

#include "J_LegL.hh"

namespace SymFunction
{

void J_LegL_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
