/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:54 GMT-05:00
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
  double t2374;
  double t2427;
  double t158;
  double t2463;
  double t2544;
  double t2548;
  double t2549;
  double t2554;
  double t2567;
  double t2569;
  double t2575;
  double t2530;
  double t2535;
  double t2537;
  double t2444;
  double t2643;
  double t2644;
  double t2552;
  double t2557;
  double t2558;
  double t2576;
  double t2591;
  double t2594;
  double t2657;
  double t2662;
  double t2663;
  double t2615;
  double t2412;
  double t2446;
  double t2479;
  double t2487;
  double t2498;
  double t2665;
  double t2666;
  double t2667;
  double t2668;
  double t2670;
  double t2671;
  double t2672;
  double t2674;
  double t2677;
  double t2646;
  double t2648;
  double t2651;
  double t2653;
  double t2654;
  double t2691;
  double t2696;
  double t2697;
  double t2600;
  double t2703;
  double t2688;
  double t2689;
  double t2698;
  double t2699;
  double t2700;
  double t2701;
  double t2702;
  double t2704;
  double t2705;
  double t2706;
  double t2707;
  double t2649;
  double t2655;
  double t2678;
  double t2710;
  double t2711;
  double t2712;
  double t2715;
  double t2719;
  double t2724;
  double t2726;
  double t2728;
  double t2734;
  double t2735;
  double t2737;
  double t2738;
  double t2740;
  double t2741;
  double t2742;
  double t2744;
  double t2745;
  double t2746;
  double t2607;
  double t2754;
  double t2757;
  double t2766;
  double t2767;
  double t2771;
  double t2788;
  double t2789;
  double t2790;
  double t2792;
  double t2793;
  double t2794;
  double t2799;
  double t2762;
  double t2774;
  double t2775;
  double t2779;
  double t2782;
  double t2783;
  double t2785;
  double t2683;
  double t2684;
  double t2685;
  double t2462;
  double t2528;
  double t2561;
  double t2595;
  double t2612;
  double t2613;
  double t2616;
  double t2622;
  double t2634;
  double t2637;
  double t2821;
  double t2851;
  double t2852;
  double t2853;
  double t2855;
  double t2868;
  double t2869;
  double t2872;
  double t2874;
  double t2865;
  double t2866;
  double t2867;
  double t2878;
  double t2879;
  double t2880;
  double t2854;
  double t2856;
  double t2857;
  double t2860;
  double t2861;
  double t2863;
  double t2873;
  double t2875;
  double t2876;
  double t2881;
  double t2882;
  double t2883;
  double t2886;
  double t2903;
  double t2904;
  double t2905;
  double t2859;
  double t2864;
  double t2877;
  double t2884;
  double t2889;
  double t2890;
  double t2891;
  double t2893;
  double t2896;
  double t2897;
  double t2898;
  double t2899;
  double t2908;
  double t2909;
  double t2919;
  double t2920;
  double t2921;
  double t2923;
  double t2924;
  double t2901;
  double t2902;
  double t2906;
  double t2907;
  double t2910;
  double t2911;
  double t2912;
  double t2913;
  double t2914;
  double t2915;
  double t2944;
  double t2945;
  double t2946;
  double t2937;
  double t2932;
  double t2933;
  double t2934;
  double t2938;
  double t2939;
  double t2940;
  double t2941;
  double t2948;
  double t2949;
  double t2922;
  double t2925;
  double t2926;
  double t2927;
  double t2928;
  double t2929;
  double t2930;
  double t2900;
  double t2916;
  double t2917;
  double t2957;
  double t2958;
  double t2959;
  double t2960;
  double t2961;
  double t2962;
  double t2963;
  double t2964;
  double t2965;
  double t2972;
  double t2979;
  double t2982;
  double t2983;
  double t2985;
  double t2986;
  double t2973;
  double t2976;
  double t2977;
  double t2918;
  double t2984;
  double t2987;
  double t2988;
  double t2952;
  double t2990;
  double t2991;
  double t2992;
  double t3001;
  double t3002;
  double t3003;
  double t3004;
  double t3005;
  double t3006;
  double t3007;
  double t2682;
  double t2708;
  double t2709;
  double t2786;
  double t2800;
  double t2801;
  double t3041;
  double t3043;
  double t3046;
  double t2818;
  double t2823;
  double t2824;
  double t2827;
  double t2830;
  double t2831;
  double t2833;
  double t2834;
  double t2836;
  double t3055;
  double t3056;
  double t3057;
  double t3058;
  double t3059;
  double t3060;
  double t3061;
  double t3062;
  double t3063;
  double t3064;
  double t3065;
  double t2805;
  double t2807;
  double t2808;
  double t2809;
  double t2810;
  double t2811;
  double t2813;
  double t2814;
  double t2815;
  double t2817;
  double t2820;
  double t3079;
  double t3080;
  double t3081;
  double t2942;
  double t2943;
  double t2947;
  double t2950;
  double t2951;
  double t2953;
  double t2954;
  double t2955;
  double t2956;
  double t2969;
  double t2970;
  double t2971;
  double t2989;
  double t2993;
  double t2994;
  double t2999;
  double t3108;
  double t3109;
  double t3110;
  double t2998;
  double t3000;
  double t3009;
  double t3010;
  double t3011;
  double t3013;
  double t3014;
  double t3130;
  double t3131;
  double t3132;
  double t3099;
  double t3100;
  double t3101;
  double t2931;
  double t2966;
  double t2967;
  t2374 = Cos(var1[3]);
  t2427 = Sin(var1[3]);
  t158 = Sin(var1[2]);
  t2463 = Cos(var1[2]);
  t2544 = Cos(var1[4]);
  t2548 = -1.*t2544;
  t2549 = 1. + t2548;
  t2554 = Sin(var1[4]);
  t2567 = -1.*t2463*t2374;
  t2569 = -1.*t158*t2427;
  t2575 = t2567 + t2569;
  t2530 = t2374*t158;
  t2535 = -1.*t2463*t2427;
  t2537 = t2530 + t2535;
  t2444 = -0.0695*t2427;
  t2643 = -1.*t2374;
  t2644 = 1. + t2643;
  t2552 = -0.0265*t2549;
  t2557 = -0.2375*t2554;
  t2558 = t2552 + t2557;
  t2576 = -0.2375*t2549;
  t2591 = 0.0265*t2554;
  t2594 = t2576 + t2591;
  t2657 = t2463*t2374;
  t2662 = t158*t2427;
  t2663 = t2657 + t2662;
  t2615 = t2544*t2537;
  t2412 = 0.0265*t2374;
  t2446 = t2412 + t2444;
  t2479 = -0.0695*t2374;
  t2487 = -0.0265*t2427;
  t2498 = t2479 + t2487;
  t2665 = -1.*t2663*t2558;
  t2666 = -1.*t2537*t2594;
  t2667 = t2544*t2663;
  t2668 = t2537*t2554;
  t2670 = t2667 + t2668;
  t2671 = 0.0265*t2670;
  t2672 = -1.*t2663*t2554;
  t2674 = t2615 + t2672;
  t2677 = -0.0325*t2674;
  t2646 = -0.0265*t2644;
  t2648 = t2646 + t2444;
  t2651 = -0.0695*t2644;
  t2653 = 0.0265*t2427;
  t2654 = t2651 + t2653;
  t2691 = -1.*t2374*t158;
  t2696 = t2463*t2427;
  t2697 = t2691 + t2696;
  t2600 = t2544*t2575;
  t2703 = t2544*t2697;
  t2688 = t158*t2648;
  t2689 = -1.*t2463*t2654;
  t2698 = -1.*t2697*t2558;
  t2699 = -1.*t2663*t2594;
  t2700 = -1.*t2697*t2554;
  t2701 = t2667 + t2700;
  t2702 = -0.0325*t2701;
  t2704 = t2663*t2554;
  t2705 = t2703 + t2704;
  t2706 = 0.0265*t2705;
  t2707 = t2688 + t2689 + t2698 + t2699 + t2702 + t2706;
  t2649 = -1.*t2463*t2648;
  t2655 = -1.*t158*t2654;
  t2678 = t2649 + t2655 + t2665 + t2666 + t2671 + t2677;
  t2710 = t2463*t2648;
  t2711 = t158*t2654;
  t2712 = -1.*t2575*t2558;
  t2715 = -1.*t2697*t2594;
  t2719 = t2697*t2554;
  t2724 = t2600 + t2719;
  t2726 = 0.0265*t2724;
  t2728 = -1.*t2575*t2554;
  t2734 = t2703 + t2728;
  t2735 = -0.0325*t2734;
  t2737 = t2710 + t2711 + t2712 + t2715 + t2726 + t2735;
  t2738 = 2.*t2737*t2707;
  t2740 = 2.*t2678*t2707;
  t2741 = t2738 + t2740;
  t2742 = Power(t2678,2);
  t2744 = Power(t2707,2);
  t2745 = 0.00085849 + t2742 + t2744;
  t2746 = Power(t2745,-1.5);
  t2607 = -1.*t2537*t2554;
  t2754 = 0.0265*t2544;
  t2757 = t2754 + t2557;
  t2766 = -0.2375*t2544;
  t2767 = -0.0265*t2554;
  t2771 = t2766 + t2767;
  t2788 = -1.*t2663*t2757;
  t2789 = -1.*t2697*t2771;
  t2790 = 0.0265*t2701;
  t2792 = -1.*t2544*t2697;
  t2793 = t2792 + t2672;
  t2794 = -0.0325*t2793;
  t2799 = t2788 + t2789 + t2790 + t2794;
  t2762 = -1.*t2537*t2757;
  t2774 = -1.*t2663*t2771;
  t2775 = -1.*t2544*t2663;
  t2779 = t2775 + t2607;
  t2782 = -0.0325*t2779;
  t2783 = 0.0265*t2674;
  t2785 = t2762 + t2774 + t2782 + t2783;
  t2683 = -1.*t2463*t2446;
  t2684 = t158*t2498;
  t2685 = t2683 + t2684 + t2665 + t2666 + t2671 + t2677;
  t2462 = -1.*t158*t2446;
  t2528 = -1.*t2463*t2498;
  t2561 = -1.*t2537*t2558;
  t2595 = -1.*t2575*t2594;
  t2612 = t2600 + t2607;
  t2613 = -0.0325*t2612;
  t2616 = t2575*t2554;
  t2622 = t2615 + t2616;
  t2634 = 0.0265*t2622;
  t2637 = t2462 + t2528 + t2561 + t2595 + t2613 + t2634;
  t2821 = 1/Sqrt(t2745);
  t2851 = Cos(var1[5]);
  t2852 = -1.*t2851;
  t2853 = 1. + t2852;
  t2855 = Sin(var1[5]);
  t2868 = Cos(var1[6]);
  t2869 = -1.*t2868;
  t2872 = 1. + t2869;
  t2874 = Sin(var1[6]);
  t2865 = t2463*t2851;
  t2866 = -1.*t158*t2855;
  t2867 = t2865 + t2866;
  t2878 = -1.*t2851*t158;
  t2879 = -1.*t2463*t2855;
  t2880 = t2878 + t2879;
  t2854 = -0.0695*t2853;
  t2856 = -0.0265*t2855;
  t2857 = t2854 + t2856;
  t2860 = -0.0265*t2853;
  t2861 = 0.0695*t2855;
  t2863 = t2860 + t2861;
  t2873 = -0.2375*t2872;
  t2875 = -0.0265*t2874;
  t2876 = t2873 + t2875;
  t2881 = -0.0265*t2872;
  t2882 = 0.2375*t2874;
  t2883 = t2881 + t2882;
  t2886 = t2868*t2867;
  t2903 = t2851*t158;
  t2904 = t2463*t2855;
  t2905 = t2903 + t2904;
  t2859 = -1.*t2463*t2857;
  t2864 = t158*t2863;
  t2877 = -1.*t2867*t2876;
  t2884 = -1.*t2880*t2883;
  t2889 = t2880*t2874;
  t2890 = t2886 + t2889;
  t2891 = -0.0325*t2890;
  t2893 = t2868*t2880;
  t2896 = -1.*t2867*t2874;
  t2897 = t2893 + t2896;
  t2898 = 0.0265*t2897;
  t2899 = t2859 + t2864 + t2877 + t2884 + t2891 + t2898;
  t2908 = -1.*t2905*t2874;
  t2909 = t2886 + t2908;
  t2919 = -0.0265*t2868;
  t2920 = -0.2375*t2874;
  t2921 = t2919 + t2920;
  t2923 = 0.2375*t2868;
  t2924 = t2923 + t2875;
  t2901 = -1.*t158*t2857;
  t2902 = -1.*t2463*t2863;
  t2906 = -1.*t2905*t2876;
  t2907 = -1.*t2867*t2883;
  t2910 = 0.0265*t2909;
  t2911 = t2868*t2905;
  t2912 = t2867*t2874;
  t2913 = t2911 + t2912;
  t2914 = -0.0325*t2913;
  t2915 = t2901 + t2902 + t2906 + t2907 + t2910 + t2914;
  t2944 = -1.*t2463*t2851;
  t2945 = t158*t2855;
  t2946 = t2944 + t2945;
  t2937 = -1.*t2880*t2874;
  t2932 = -1.*t2867*t2921;
  t2933 = -1.*t2880*t2924;
  t2934 = -1.*t2868*t2867;
  t2938 = t2934 + t2937;
  t2939 = 0.0265*t2938;
  t2940 = -0.0325*t2897;
  t2941 = t2932 + t2933 + t2939 + t2940;
  t2948 = t2868*t2946;
  t2949 = t2948 + t2937;
  t2922 = -1.*t2905*t2921;
  t2925 = -1.*t2867*t2924;
  t2926 = -0.0325*t2909;
  t2927 = -1.*t2868*t2905;
  t2928 = t2927 + t2896;
  t2929 = 0.0265*t2928;
  t2930 = t2922 + t2925 + t2926 + t2929;
  t2900 = Power(t2899,2);
  t2916 = Power(t2915,2);
  t2917 = 0.00085849 + t2900 + t2916;
  t2957 = t158*t2857;
  t2958 = t2463*t2863;
  t2959 = -1.*t2880*t2876;
  t2960 = -1.*t2946*t2883;
  t2961 = 0.0265*t2949;
  t2962 = t2946*t2874;
  t2963 = t2893 + t2962;
  t2964 = -0.0325*t2963;
  t2965 = t2957 + t2958 + t2959 + t2960 + t2961 + t2964;
  t2972 = Power(t2917,-1.5);
  t2979 = -0.0265*t2851;
  t2982 = -0.0695*t2855;
  t2983 = t2979 + t2982;
  t2985 = 0.0695*t2851;
  t2986 = t2985 + t2856;
  t2973 = 2.*t2899*t2915;
  t2976 = 2.*t2899*t2965;
  t2977 = t2973 + t2976;
  t2918 = 1/Sqrt(t2917);
  t2984 = -1.*t158*t2983;
  t2987 = -1.*t2463*t2986;
  t2988 = t2984 + t2987 + t2877 + t2884 + t2891 + t2898;
  t2952 = -1.*t2946*t2874;
  t2990 = -1.*t2463*t2983;
  t2991 = t158*t2986;
  t2992 = t2990 + t2991 + t2959 + t2960 + t2961 + t2964;
  t3001 = -1.*t2946*t2876;
  t3002 = -1.*t2905*t2883;
  t3003 = t2905*t2874;
  t3004 = t2948 + t3003;
  t3005 = -0.0325*t3004;
  t3006 = t2911 + t2952;
  t3007 = 0.0265*t3006;
  t2682 = 2.*t2637*t2678;
  t2708 = 2.*t2685*t2707;
  t2709 = t2682 + t2708;
  t2786 = 2.*t2678*t2785;
  t2800 = 2.*t2799*t2707;
  t2801 = t2786 + t2800;
  t3041 = -0.0265*t2374;
  t3043 = 0.0695*t2427;
  t3046 = t3041 + t3043;
  t2818 = 2.*t2785*t2707;
  t2823 = 2.*t2737*t2685;
  t2824 = 2.*t2685*t2678;
  t2827 = 2.*t2637*t2707;
  t2830 = t158*t2446;
  t2831 = t2463*t2498;
  t2833 = t2830 + t2831 + t2698 + t2699 + t2702 + t2706;
  t2834 = 2.*t2833*t2707;
  t2836 = t2823 + t2824 + t2827 + t2834;
  t3055 = -1.*t2575*t2757;
  t3056 = -1.*t2537*t2771;
  t3057 = 0.0265*t2612;
  t3058 = -1.*t2544*t2537;
  t3059 = t3058 + t2728;
  t3060 = -0.0325*t3059;
  t3061 = t3055 + t3056 + t3057 + t3060;
  t3062 = 2.*t3061*t2678;
  t3063 = 2.*t2637*t2785;
  t3064 = 2.*t2685*t2799;
  t3065 = t3062 + t3063 + t3064 + t2818;
  t2805 = 2.*t2737*t2799;
  t2807 = 2.*t2678*t2799;
  t2808 = -1.*t2697*t2757;
  t2809 = -1.*t2575*t2771;
  t2810 = -1.*t2544*t2575;
  t2811 = t2810 + t2700;
  t2813 = -0.0325*t2811;
  t2814 = 0.0265*t2734;
  t2815 = t2808 + t2809 + t2813 + t2814;
  t2817 = 2.*t2815*t2707;
  t2820 = t2805 + t2807 + t2817 + t2818;
  t3079 = -0.0265*t2544;
  t3080 = 0.2375*t2554;
  t3081 = t3079 + t3080;
  t2942 = 2.*t2941*t2915;
  t2943 = -1.*t2880*t2921;
  t2947 = -1.*t2946*t2924;
  t2950 = -0.0325*t2949;
  t2951 = -1.*t2868*t2880;
  t2953 = t2951 + t2952;
  t2954 = 0.0265*t2953;
  t2955 = t2943 + t2947 + t2950 + t2954;
  t2956 = 2.*t2899*t2955;
  t2969 = 2.*t2941*t2899;
  t2970 = 2.*t2930*t2915;
  t2971 = t2969 + t2970;
  t2989 = 2.*t2988*t2915;
  t2993 = 2.*t2899*t2992;
  t2994 = t2989 + t2993;
  t2999 = t158*t2983;
  t3108 = -0.0695*t2851;
  t3109 = 0.0265*t2855;
  t3110 = t3108 + t3109;
  t2998 = 2.*t2988*t2899;
  t3000 = t2463*t2986;
  t3009 = t2999 + t3000 + t3001 + t3002 + t3005 + t3007;
  t3010 = 2.*t2899*t3009;
  t3011 = 2.*t2915*t2992;
  t3013 = 2.*t2992*t2965;
  t3014 = t2998 + t3010 + t3011 + t3013;
  t3130 = -0.2375*t2868;
  t3131 = 0.0265*t2874;
  t3132 = t3130 + t3131;
  t3099 = 2.*t2988*t2930;
  t3100 = 2.*t2941*t2992;
  t3101 = t3099 + t2942 + t2956 + t3100;
  t2931 = 2.*t2899*t2930;
  t2966 = 2.*t2941*t2965;
  t2967 = t2931 + t2942 + t2956 + t2966;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=-0.25*Power(t2741,2)*t2746*var1[9] + 0.5*(2.*(t2561 + t2595 + t2613 + t2634 - 1.*t158*t2648 + t2463*t2654)*t2707 + 2.*t2678*t2737 + 2.*Power(t2737,2) + 2.*t2744)*t2821*var1[9] - 0.25*t2709*t2741*t2746*var1[10] + 0.5*t2821*t2836*var1[10] - 0.25*t2741*t2746*t2801*var1[11] + 0.5*t2820*t2821*var1[11];
  p_output1[5]=-0.25*t2972*Power(t2977,2)*var1[9] + 0.5*t2918*(2.*t2900 + 2.*t2915*t2965 + 2.*Power(t2965,2) + 2.*t2899*(t2463*t2857 - 1.*t158*t2863 + t3001 + t3002 + t3005 + t3007))*var1[9] - 0.25*t2972*t2977*t2994*var1[12] + 0.5*t2918*t3014*var1[12] + 0.5*t2918*t2967*var1[13] - 0.25*t2971*t2972*t2977*var1[13];
  p_output1[6]=-0.25*t2709*t2741*t2746*var1[9] + 0.5*t2821*t2836*var1[9] - 0.25*Power(t2709,2)*t2746*var1[10] + 0.5*t2821*(2.*Power(t2637,2) + 2.*Power(t2685,2) + 2.*t2707*(t2528 + t2561 + t2595 + t2613 + t2634 + t158*t3046) + 2.*t2678*(-1.*t158*t2498 + t2712 + t2715 + t2726 + t2735 - 1.*t2463*t3046))*var1[10] - 0.25*t2709*t2746*t2801*var1[11] + 0.5*t2821*t3065*var1[11];
  p_output1[7]=0;
  p_output1[8]=-0.25*t2741*t2746*t2801*var1[9] + 0.5*t2820*t2821*var1[9] - 0.25*t2709*t2746*t2801*var1[10] + 0.5*t2821*t3065*var1[10] - 0.25*t2746*Power(t2801,2)*var1[11] + 0.5*t2821*(2.*Power(t2785,2) + 2.*Power(t2799,2) + 2.*t2678*(0.0265*t2779 + t3056 - 0.0325*(t2704 + t3058) - 1.*t2663*t3081) + 2.*t2707*(t2774 - 0.0325*(t2719 + t2775) + 0.0265*t2793 - 1.*t2697*t3081))*var1[11];
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=-0.25*t2972*t2977*t2994*var1[9] + 0.5*t2918*t3014*var1[9] - 0.25*t2972*Power(t2994,2)*var1[12] + 0.5*t2918*(2.*Power(t2988,2) + 2.*Power(t2992,2) + 2.*t2915*(t2959 + t2960 + t2961 + t2964 + t2990 - 1.*t158*t3110) + 2.*t2899*(t2999 + t3001 + t3002 + t3005 + t3007 - 1.*t2463*t3110))*var1[12] - 0.25*t2971*t2972*t2994*var1[13] + 0.5*t2918*t3101*var1[13];
  p_output1[12]=0;
  p_output1[13]=0.5*t2918*t2967*var1[9] - 0.25*t2971*t2972*t2977*var1[9] - 0.25*t2971*t2972*t2994*var1[12] + 0.5*t2918*t3101*var1[12] - 0.25*Power(t2971,2)*t2972*var1[13] + 0.5*t2918*(2.*Power(t2930,2) + 2.*Power(t2941,2) + 2.*t2899*(-0.0325*t2938 + t2943 + 0.0265*(t2912 + t2951) - 1.*t2867*t3132) + 2.*t2915*(-0.0325*t2928 + t2932 + 0.0265*(t2934 + t3003) - 1.*t2905*t3132))*var1[13];
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
