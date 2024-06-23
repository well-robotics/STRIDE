/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:52 GMT-05:00
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
  double t2320;
  double t2338;
  double t2352;
  double t2371;
  double t2287;
  double t2400;
  double t2436;
  double t2438;
  double t2440;
  double t2442;
  double t2429;
  double t2432;
  double t2435;
  double t2448;
  double t2451;
  double t2457;
  double t2356;
  double t2374;
  double t2375;
  double t2412;
  double t2413;
  double t2426;
  double t2441;
  double t2444;
  double t2445;
  double t2462;
  double t2463;
  double t2465;
  double t2505;
  double t2507;
  double t2508;
  double t2477;
  double t2517;
  double t2544;
  double t2545;
  double t2546;
  double t2498;
  double t2500;
  double t2503;
  double t2515;
  double t2524;
  double t2528;
  double t2530;
  double t2532;
  double t2535;
  double t2537;
  double t2538;
  double t2541;
  double t2542;
  double t2543;
  double t2547;
  double t2548;
  double t2549;
  double t2551;
  double t2552;
  double t2553;
  double t2554;
  double t2555;
  double t2556;
  double t2564;
  double t2565;
  double t2566;
  double t2568;
  double t2579;
  double t2580;
  double t2581;
  double t2583;
  double t2576;
  double t2577;
  double t2578;
  double t2587;
  double t2588;
  double t2589;
  double t2567;
  double t2569;
  double t2570;
  double t2572;
  double t2573;
  double t2574;
  double t2582;
  double t2584;
  double t2585;
  double t2590;
  double t2591;
  double t2592;
  double t2594;
  double t2606;
  double t2607;
  double t2608;
  double t2571;
  double t2575;
  double t2586;
  double t2593;
  double t2595;
  double t2596;
  double t2597;
  double t2598;
  double t2599;
  double t2600;
  double t2601;
  double t2602;
  double t2604;
  double t2605;
  double t2609;
  double t2610;
  double t2611;
  double t2612;
  double t2613;
  double t2615;
  double t2616;
  double t2617;
  double t2620;
  double t2621;
  double t2633;
  double t2634;
  double t2635;
  double t2470;
  double t2651;
  double t2652;
  double t2654;
  double t2655;
  double t2656;
  double t2559;
  double t2560;
  double t2561;
  double t2562;
  double t2664;
  double t2685;
  double t2686;
  double t2688;
  double t2689;
  double t2690;
  double t2603;
  double t2622;
  double t2623;
  double t2625;
  double t2711;
  double t2712;
  double t2713;
  double t2716;
  double t2718;
  double t2630;
  double t2636;
  double t2637;
  double t2638;
  double t2641;
  double t2642;
  double t2643;
  double t2644;
  double t2645;
  double t2728;
  double t2729;
  double t2730;
  double t2732;
  double t2733;
  t2320 = Cos(var1[3]);
  t2338 = -1.*t2320;
  t2352 = 1. + t2338;
  t2371 = Sin(var1[3]);
  t2287 = Cos(var1[2]);
  t2400 = Sin(var1[2]);
  t2436 = Cos(var1[4]);
  t2438 = -1.*t2436;
  t2440 = 1. + t2438;
  t2442 = Sin(var1[4]);
  t2429 = -1.*t2287*t2320;
  t2432 = -1.*t2400*t2371;
  t2435 = t2429 + t2432;
  t2448 = -1.*t2320*t2400;
  t2451 = t2287*t2371;
  t2457 = t2448 + t2451;
  t2356 = -0.0265*t2352;
  t2374 = -0.0695*t2371;
  t2375 = t2356 + t2374;
  t2412 = -0.0695*t2352;
  t2413 = 0.0265*t2371;
  t2426 = t2412 + t2413;
  t2441 = -0.0265*t2440;
  t2444 = -0.2375*t2442;
  t2445 = t2441 + t2444;
  t2462 = -0.2375*t2440;
  t2463 = 0.0265*t2442;
  t2465 = t2462 + t2463;
  t2505 = t2287*t2320;
  t2507 = t2400*t2371;
  t2508 = t2505 + t2507;
  t2477 = t2436*t2457;
  t2517 = t2436*t2508;
  t2544 = t2320*t2400;
  t2545 = -1.*t2287*t2371;
  t2546 = t2544 + t2545;
  t2498 = t2400*t2375;
  t2500 = -1.*t2287*t2426;
  t2503 = -1.*t2457*t2445;
  t2515 = -1.*t2508*t2465;
  t2524 = -1.*t2457*t2442;
  t2528 = t2517 + t2524;
  t2530 = -0.0325*t2528;
  t2532 = t2508*t2442;
  t2535 = t2477 + t2532;
  t2537 = 0.0265*t2535;
  t2538 = t2498 + t2500 + t2503 + t2515 + t2530 + t2537;
  t2541 = -1.*t2287*t2375;
  t2542 = -1.*t2400*t2426;
  t2543 = -1.*t2508*t2445;
  t2547 = -1.*t2546*t2465;
  t2548 = t2546*t2442;
  t2549 = t2517 + t2548;
  t2551 = 0.0265*t2549;
  t2552 = t2436*t2546;
  t2553 = -1.*t2508*t2442;
  t2554 = t2552 + t2553;
  t2555 = -0.0325*t2554;
  t2556 = t2541 + t2542 + t2543 + t2547 + t2551 + t2555;
  t2564 = Cos(var1[5]);
  t2565 = -1.*t2564;
  t2566 = 1. + t2565;
  t2568 = Sin(var1[5]);
  t2579 = Cos(var1[6]);
  t2580 = -1.*t2579;
  t2581 = 1. + t2580;
  t2583 = Sin(var1[6]);
  t2576 = t2287*t2564;
  t2577 = -1.*t2400*t2568;
  t2578 = t2576 + t2577;
  t2587 = -1.*t2564*t2400;
  t2588 = -1.*t2287*t2568;
  t2589 = t2587 + t2588;
  t2567 = -0.0695*t2566;
  t2569 = -0.0265*t2568;
  t2570 = t2567 + t2569;
  t2572 = -0.0265*t2566;
  t2573 = 0.0695*t2568;
  t2574 = t2572 + t2573;
  t2582 = -0.2375*t2581;
  t2584 = -0.0265*t2583;
  t2585 = t2582 + t2584;
  t2590 = -0.0265*t2581;
  t2591 = 0.2375*t2583;
  t2592 = t2590 + t2591;
  t2594 = t2579*t2578;
  t2606 = t2564*t2400;
  t2607 = t2287*t2568;
  t2608 = t2606 + t2607;
  t2571 = -1.*t2287*t2570;
  t2575 = t2400*t2574;
  t2586 = -1.*t2578*t2585;
  t2593 = -1.*t2589*t2592;
  t2595 = t2589*t2583;
  t2596 = t2594 + t2595;
  t2597 = -0.0325*t2596;
  t2598 = t2579*t2589;
  t2599 = -1.*t2578*t2583;
  t2600 = t2598 + t2599;
  t2601 = 0.0265*t2600;
  t2602 = t2571 + t2575 + t2586 + t2593 + t2597 + t2601;
  t2604 = -1.*t2400*t2570;
  t2605 = -1.*t2287*t2574;
  t2609 = -1.*t2608*t2585;
  t2610 = -1.*t2578*t2592;
  t2611 = -1.*t2608*t2583;
  t2612 = t2594 + t2611;
  t2613 = 0.0265*t2612;
  t2615 = t2579*t2608;
  t2616 = t2578*t2583;
  t2617 = t2615 + t2616;
  t2620 = -0.0325*t2617;
  t2621 = t2604 + t2605 + t2609 + t2610 + t2613 + t2620;
  t2633 = -1.*t2287*t2564;
  t2634 = t2400*t2568;
  t2635 = t2633 + t2634;
  t2470 = t2436*t2435;
  t2651 = 0.0265*t2320;
  t2652 = t2651 + t2374;
  t2654 = -0.0695*t2320;
  t2655 = -0.0265*t2371;
  t2656 = t2654 + t2655;
  t2559 = Power(t2556,2);
  t2560 = Power(t2538,2);
  t2561 = 0.00085849 + t2559 + t2560;
  t2562 = 1/Sqrt(t2561);
  t2664 = -1.*t2546*t2442;
  t2685 = 0.0265*t2436;
  t2686 = t2685 + t2444;
  t2688 = -0.2375*t2436;
  t2689 = -0.0265*t2442;
  t2690 = t2688 + t2689;
  t2603 = Power(t2602,2);
  t2622 = Power(t2621,2);
  t2623 = 0.00085849 + t2603 + t2622;
  t2625 = 1/Sqrt(t2623);
  t2711 = -0.0265*t2564;
  t2712 = -0.0695*t2568;
  t2713 = t2711 + t2712;
  t2716 = 0.0695*t2564;
  t2718 = t2716 + t2569;
  t2630 = -1.*t2589*t2585;
  t2636 = -1.*t2635*t2592;
  t2637 = t2579*t2635;
  t2638 = -1.*t2589*t2583;
  t2641 = t2637 + t2638;
  t2642 = 0.0265*t2641;
  t2643 = t2635*t2583;
  t2644 = t2598 + t2643;
  t2645 = -0.0325*t2644;
  t2728 = -0.0265*t2579;
  t2729 = -0.2375*t2583;
  t2730 = t2728 + t2729;
  t2732 = 0.2375*t2579;
  t2733 = t2732 + t2584;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0.5*(2.*(t2287*t2375 + t2400*t2426 - 1.*t2435*t2445 - 1.*t2457*t2465 + 0.0265*(t2442*t2457 + t2470) - 0.0325*(-1.*t2435*t2442 + t2477))*t2538 + 2.*t2538*t2556)*t2562;
  p_output1[5]=0.5*t2625*(2.*t2602*t2621 + 2.*t2602*(t2400*t2570 + t2287*t2574 + t2630 + t2636 + t2642 + t2645));
  p_output1[6]=0.5*t2562*(2.*t2538*(t2543 + t2547 + t2551 + t2555 - 1.*t2287*t2652 + t2400*t2656) + 2.*t2556*(-1.*t2435*t2465 - 1.*t2445*t2546 + 0.0265*(t2435*t2442 + t2552) - 1.*t2400*t2652 - 1.*t2287*t2656 - 0.0325*(t2470 + t2664)));
  p_output1[7]=0;
  p_output1[8]=0.5*t2562*(2.*t2538*(0.0265*t2528 - 0.0325*(-1.*t2436*t2457 + t2553) - 1.*t2508*t2686 - 1.*t2457*t2690) + 2.*t2556*(0.0265*t2554 - 0.0325*(-1.*t2436*t2508 + t2664) - 1.*t2546*t2686 - 1.*t2508*t2690));
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=0.5*t2625*(2.*t2621*(t2586 + t2593 + t2597 + t2601 - 1.*t2400*t2713 - 1.*t2287*t2718) + 2.*t2602*(t2630 + t2636 + t2642 + t2645 - 1.*t2287*t2713 + t2400*t2718));
  p_output1[12]=0;
  p_output1[13]=0.5*t2625*(2.*t2621*(0.0265*(t2599 - 1.*t2579*t2608) - 0.0325*t2612 - 1.*t2608*t2730 - 1.*t2578*t2733) + 2.*t2602*(-0.0325*t2600 + 0.0265*(-1.*t2578*t2579 + t2638) - 1.*t2578*t2730 - 1.*t2589*t2733));
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
