/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:08 GMT-05:00
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
  double t2445;
  double t2497;
  double t744;
  double t2522;
  double t2604;
  double t2605;
  double t2608;
  double t2613;
  double t2625;
  double t2628;
  double t2632;
  double t2588;
  double t2593;
  double t2595;
  double t2502;
  double t2700;
  double t2702;
  double t2610;
  double t2614;
  double t2615;
  double t2638;
  double t2650;
  double t2651;
  double t2714;
  double t2719;
  double t2721;
  double t2672;
  double t2469;
  double t2507;
  double t2539;
  double t2551;
  double t2561;
  double t2722;
  double t2723;
  double t2724;
  double t2725;
  double t2727;
  double t2728;
  double t2729;
  double t2731;
  double t2734;
  double t2703;
  double t2705;
  double t2709;
  double t2710;
  double t2711;
  double t2748;
  double t2753;
  double t2754;
  double t2659;
  double t2760;
  double t2745;
  double t2747;
  double t2755;
  double t2756;
  double t2757;
  double t2758;
  double t2759;
  double t2761;
  double t2762;
  double t2763;
  double t2764;
  double t2706;
  double t2713;
  double t2735;
  double t2767;
  double t2768;
  double t2770;
  double t2772;
  double t2776;
  double t2781;
  double t2783;
  double t2785;
  double t2791;
  double t2792;
  double t2794;
  double t2795;
  double t2797;
  double t2798;
  double t2799;
  double t2801;
  double t2802;
  double t2803;
  double t2668;
  double t2811;
  double t2814;
  double t2823;
  double t2824;
  double t2828;
  double t2845;
  double t2846;
  double t2847;
  double t2849;
  double t2850;
  double t2851;
  double t2856;
  double t2819;
  double t2831;
  double t2832;
  double t2836;
  double t2839;
  double t2840;
  double t2842;
  double t2740;
  double t2741;
  double t2743;
  double t2519;
  double t2585;
  double t2618;
  double t2652;
  double t2669;
  double t2670;
  double t2673;
  double t2679;
  double t2693;
  double t2695;
  double t2878;
  double t2908;
  double t2909;
  double t2910;
  double t2912;
  double t2925;
  double t2926;
  double t2929;
  double t2931;
  double t2922;
  double t2923;
  double t2924;
  double t2935;
  double t2936;
  double t2937;
  double t2911;
  double t2913;
  double t2914;
  double t2917;
  double t2918;
  double t2920;
  double t2930;
  double t2932;
  double t2933;
  double t2938;
  double t2939;
  double t2940;
  double t2943;
  double t2960;
  double t2961;
  double t2962;
  double t2916;
  double t2921;
  double t2934;
  double t2941;
  double t2946;
  double t2947;
  double t2948;
  double t2950;
  double t2953;
  double t2954;
  double t2955;
  double t2956;
  double t2965;
  double t2966;
  double t2976;
  double t2977;
  double t2978;
  double t2980;
  double t2981;
  double t2958;
  double t2959;
  double t2963;
  double t2964;
  double t2967;
  double t2968;
  double t2969;
  double t2970;
  double t2971;
  double t2972;
  double t3001;
  double t3002;
  double t3003;
  double t2994;
  double t2989;
  double t2990;
  double t2991;
  double t2995;
  double t2996;
  double t2997;
  double t2998;
  double t3005;
  double t3006;
  double t2979;
  double t2982;
  double t2983;
  double t2984;
  double t2985;
  double t2986;
  double t2987;
  double t2957;
  double t2973;
  double t2974;
  double t3014;
  double t3015;
  double t3016;
  double t3017;
  double t3018;
  double t3019;
  double t3020;
  double t3021;
  double t3022;
  double t3029;
  double t3036;
  double t3039;
  double t3040;
  double t3042;
  double t3043;
  double t3030;
  double t3033;
  double t3034;
  double t2975;
  double t3041;
  double t3044;
  double t3045;
  double t3009;
  double t3047;
  double t3048;
  double t3049;
  double t3058;
  double t3059;
  double t3060;
  double t3061;
  double t3062;
  double t3063;
  double t3064;
  double t2739;
  double t2765;
  double t2766;
  double t2843;
  double t2857;
  double t2858;
  double t3098;
  double t3100;
  double t3103;
  double t2875;
  double t2880;
  double t2881;
  double t2884;
  double t2887;
  double t2888;
  double t2890;
  double t2891;
  double t2893;
  double t3112;
  double t3113;
  double t3114;
  double t3115;
  double t3116;
  double t3117;
  double t3118;
  double t3119;
  double t3120;
  double t3121;
  double t3122;
  double t2862;
  double t2864;
  double t2865;
  double t2866;
  double t2867;
  double t2868;
  double t2870;
  double t2871;
  double t2872;
  double t2874;
  double t2877;
  double t3136;
  double t3137;
  double t3138;
  double t2999;
  double t3000;
  double t3004;
  double t3007;
  double t3008;
  double t3010;
  double t3011;
  double t3012;
  double t3013;
  double t3026;
  double t3027;
  double t3028;
  double t3046;
  double t3050;
  double t3051;
  double t3056;
  double t3165;
  double t3166;
  double t3167;
  double t3055;
  double t3057;
  double t3066;
  double t3067;
  double t3068;
  double t3070;
  double t3071;
  double t3187;
  double t3188;
  double t3189;
  double t3156;
  double t3157;
  double t3158;
  double t2988;
  double t3023;
  double t3024;
  t2445 = Cos(var1[3]);
  t2497 = Sin(var1[3]);
  t744 = Sin(var1[2]);
  t2522 = Cos(var1[2]);
  t2604 = Cos(var1[4]);
  t2605 = -1.*t2604;
  t2608 = 1. + t2605;
  t2613 = Sin(var1[4]);
  t2625 = -1.*t2522*t2445;
  t2628 = -1.*t744*t2497;
  t2632 = t2625 + t2628;
  t2588 = t2445*t744;
  t2593 = -1.*t2522*t2497;
  t2595 = t2588 + t2593;
  t2502 = -0.0695*t2497;
  t2700 = -1.*t2445;
  t2702 = 1. + t2700;
  t2610 = -0.0265*t2608;
  t2614 = -0.2375*t2613;
  t2615 = t2610 + t2614;
  t2638 = -0.2375*t2608;
  t2650 = 0.0265*t2613;
  t2651 = t2638 + t2650;
  t2714 = t2522*t2445;
  t2719 = t744*t2497;
  t2721 = t2714 + t2719;
  t2672 = t2604*t2595;
  t2469 = 0.0265*t2445;
  t2507 = t2469 + t2502;
  t2539 = -0.0695*t2445;
  t2551 = -0.0265*t2497;
  t2561 = t2539 + t2551;
  t2722 = -1.*t2721*t2615;
  t2723 = -1.*t2595*t2651;
  t2724 = t2604*t2721;
  t2725 = t2595*t2613;
  t2727 = t2724 + t2725;
  t2728 = 0.0265*t2727;
  t2729 = -1.*t2721*t2613;
  t2731 = t2672 + t2729;
  t2734 = 0.0115*t2731;
  t2703 = -0.0265*t2702;
  t2705 = t2703 + t2502;
  t2709 = -0.0695*t2702;
  t2710 = 0.0265*t2497;
  t2711 = t2709 + t2710;
  t2748 = -1.*t2445*t744;
  t2753 = t2522*t2497;
  t2754 = t2748 + t2753;
  t2659 = t2604*t2632;
  t2760 = t2604*t2754;
  t2745 = t744*t2705;
  t2747 = -1.*t2522*t2711;
  t2755 = -1.*t2754*t2615;
  t2756 = -1.*t2721*t2651;
  t2757 = -1.*t2754*t2613;
  t2758 = t2724 + t2757;
  t2759 = 0.0115*t2758;
  t2761 = t2721*t2613;
  t2762 = t2760 + t2761;
  t2763 = 0.0265*t2762;
  t2764 = t2745 + t2747 + t2755 + t2756 + t2759 + t2763;
  t2706 = -1.*t2522*t2705;
  t2713 = -1.*t744*t2711;
  t2735 = t2706 + t2713 + t2722 + t2723 + t2728 + t2734;
  t2767 = t2522*t2705;
  t2768 = t744*t2711;
  t2770 = -1.*t2632*t2615;
  t2772 = -1.*t2754*t2651;
  t2776 = t2754*t2613;
  t2781 = t2659 + t2776;
  t2783 = 0.0265*t2781;
  t2785 = -1.*t2632*t2613;
  t2791 = t2760 + t2785;
  t2792 = 0.0115*t2791;
  t2794 = t2767 + t2768 + t2770 + t2772 + t2783 + t2792;
  t2795 = 2.*t2794*t2764;
  t2797 = 2.*t2735*t2764;
  t2798 = t2795 + t2797;
  t2799 = Power(t2735,2);
  t2801 = Power(t2764,2);
  t2802 = 0.00085849 + t2799 + t2801;
  t2803 = Power(t2802,-1.5);
  t2668 = -1.*t2595*t2613;
  t2811 = 0.0265*t2604;
  t2814 = t2811 + t2614;
  t2823 = -0.2375*t2604;
  t2824 = -0.0265*t2613;
  t2828 = t2823 + t2824;
  t2845 = -1.*t2721*t2814;
  t2846 = -1.*t2754*t2828;
  t2847 = 0.0265*t2758;
  t2849 = -1.*t2604*t2754;
  t2850 = t2849 + t2729;
  t2851 = 0.0115*t2850;
  t2856 = t2845 + t2846 + t2847 + t2851;
  t2819 = -1.*t2595*t2814;
  t2831 = -1.*t2721*t2828;
  t2832 = -1.*t2604*t2721;
  t2836 = t2832 + t2668;
  t2839 = 0.0115*t2836;
  t2840 = 0.0265*t2731;
  t2842 = t2819 + t2831 + t2839 + t2840;
  t2740 = -1.*t2522*t2507;
  t2741 = t744*t2561;
  t2743 = t2740 + t2741 + t2722 + t2723 + t2728 + t2734;
  t2519 = -1.*t744*t2507;
  t2585 = -1.*t2522*t2561;
  t2618 = -1.*t2595*t2615;
  t2652 = -1.*t2632*t2651;
  t2669 = t2659 + t2668;
  t2670 = 0.0115*t2669;
  t2673 = t2632*t2613;
  t2679 = t2672 + t2673;
  t2693 = 0.0265*t2679;
  t2695 = t2519 + t2585 + t2618 + t2652 + t2670 + t2693;
  t2878 = 1/Sqrt(t2802);
  t2908 = Cos(var1[5]);
  t2909 = -1.*t2908;
  t2910 = 1. + t2909;
  t2912 = Sin(var1[5]);
  t2925 = Cos(var1[6]);
  t2926 = -1.*t2925;
  t2929 = 1. + t2926;
  t2931 = Sin(var1[6]);
  t2922 = t2522*t2908;
  t2923 = -1.*t744*t2912;
  t2924 = t2922 + t2923;
  t2935 = -1.*t2908*t744;
  t2936 = -1.*t2522*t2912;
  t2937 = t2935 + t2936;
  t2911 = -0.0695*t2910;
  t2913 = -0.0265*t2912;
  t2914 = t2911 + t2913;
  t2917 = -0.0265*t2910;
  t2918 = 0.0695*t2912;
  t2920 = t2917 + t2918;
  t2930 = -0.2375*t2929;
  t2932 = -0.0265*t2931;
  t2933 = t2930 + t2932;
  t2938 = -0.0265*t2929;
  t2939 = 0.2375*t2931;
  t2940 = t2938 + t2939;
  t2943 = t2925*t2924;
  t2960 = t2908*t744;
  t2961 = t2522*t2912;
  t2962 = t2960 + t2961;
  t2916 = -1.*t2522*t2914;
  t2921 = t744*t2920;
  t2934 = -1.*t2924*t2933;
  t2941 = -1.*t2937*t2940;
  t2946 = t2937*t2931;
  t2947 = t2943 + t2946;
  t2948 = 0.0115*t2947;
  t2950 = t2925*t2937;
  t2953 = -1.*t2924*t2931;
  t2954 = t2950 + t2953;
  t2955 = 0.0265*t2954;
  t2956 = t2916 + t2921 + t2934 + t2941 + t2948 + t2955;
  t2965 = -1.*t2962*t2931;
  t2966 = t2943 + t2965;
  t2976 = -0.0265*t2925;
  t2977 = -0.2375*t2931;
  t2978 = t2976 + t2977;
  t2980 = 0.2375*t2925;
  t2981 = t2980 + t2932;
  t2958 = -1.*t744*t2914;
  t2959 = -1.*t2522*t2920;
  t2963 = -1.*t2962*t2933;
  t2964 = -1.*t2924*t2940;
  t2967 = 0.0265*t2966;
  t2968 = t2925*t2962;
  t2969 = t2924*t2931;
  t2970 = t2968 + t2969;
  t2971 = 0.0115*t2970;
  t2972 = t2958 + t2959 + t2963 + t2964 + t2967 + t2971;
  t3001 = -1.*t2522*t2908;
  t3002 = t744*t2912;
  t3003 = t3001 + t3002;
  t2994 = -1.*t2937*t2931;
  t2989 = -1.*t2924*t2978;
  t2990 = -1.*t2937*t2981;
  t2991 = -1.*t2925*t2924;
  t2995 = t2991 + t2994;
  t2996 = 0.0265*t2995;
  t2997 = 0.0115*t2954;
  t2998 = t2989 + t2990 + t2996 + t2997;
  t3005 = t2925*t3003;
  t3006 = t3005 + t2994;
  t2979 = -1.*t2962*t2978;
  t2982 = -1.*t2924*t2981;
  t2983 = 0.0115*t2966;
  t2984 = -1.*t2925*t2962;
  t2985 = t2984 + t2953;
  t2986 = 0.0265*t2985;
  t2987 = t2979 + t2982 + t2983 + t2986;
  t2957 = Power(t2956,2);
  t2973 = Power(t2972,2);
  t2974 = 0.00085849 + t2957 + t2973;
  t3014 = t744*t2914;
  t3015 = t2522*t2920;
  t3016 = -1.*t2937*t2933;
  t3017 = -1.*t3003*t2940;
  t3018 = 0.0265*t3006;
  t3019 = t3003*t2931;
  t3020 = t2950 + t3019;
  t3021 = 0.0115*t3020;
  t3022 = t3014 + t3015 + t3016 + t3017 + t3018 + t3021;
  t3029 = Power(t2974,-1.5);
  t3036 = -0.0265*t2908;
  t3039 = -0.0695*t2912;
  t3040 = t3036 + t3039;
  t3042 = 0.0695*t2908;
  t3043 = t3042 + t2913;
  t3030 = 2.*t2956*t2972;
  t3033 = 2.*t2956*t3022;
  t3034 = t3030 + t3033;
  t2975 = 1/Sqrt(t2974);
  t3041 = -1.*t744*t3040;
  t3044 = -1.*t2522*t3043;
  t3045 = t3041 + t3044 + t2934 + t2941 + t2948 + t2955;
  t3009 = -1.*t3003*t2931;
  t3047 = -1.*t2522*t3040;
  t3048 = t744*t3043;
  t3049 = t3047 + t3048 + t3016 + t3017 + t3018 + t3021;
  t3058 = -1.*t3003*t2933;
  t3059 = -1.*t2962*t2940;
  t3060 = t2962*t2931;
  t3061 = t3005 + t3060;
  t3062 = 0.0115*t3061;
  t3063 = t2968 + t3009;
  t3064 = 0.0265*t3063;
  t2739 = 2.*t2695*t2735;
  t2765 = 2.*t2743*t2764;
  t2766 = t2739 + t2765;
  t2843 = 2.*t2735*t2842;
  t2857 = 2.*t2856*t2764;
  t2858 = t2843 + t2857;
  t3098 = -0.0265*t2445;
  t3100 = 0.0695*t2497;
  t3103 = t3098 + t3100;
  t2875 = 2.*t2842*t2764;
  t2880 = 2.*t2794*t2743;
  t2881 = 2.*t2743*t2735;
  t2884 = 2.*t2695*t2764;
  t2887 = t744*t2507;
  t2888 = t2522*t2561;
  t2890 = t2887 + t2888 + t2755 + t2756 + t2759 + t2763;
  t2891 = 2.*t2890*t2764;
  t2893 = t2880 + t2881 + t2884 + t2891;
  t3112 = -1.*t2632*t2814;
  t3113 = -1.*t2595*t2828;
  t3114 = 0.0265*t2669;
  t3115 = -1.*t2604*t2595;
  t3116 = t3115 + t2785;
  t3117 = 0.0115*t3116;
  t3118 = t3112 + t3113 + t3114 + t3117;
  t3119 = 2.*t3118*t2735;
  t3120 = 2.*t2695*t2842;
  t3121 = 2.*t2743*t2856;
  t3122 = t3119 + t3120 + t3121 + t2875;
  t2862 = 2.*t2794*t2856;
  t2864 = 2.*t2735*t2856;
  t2865 = -1.*t2754*t2814;
  t2866 = -1.*t2632*t2828;
  t2867 = -1.*t2604*t2632;
  t2868 = t2867 + t2757;
  t2870 = 0.0115*t2868;
  t2871 = 0.0265*t2791;
  t2872 = t2865 + t2866 + t2870 + t2871;
  t2874 = 2.*t2872*t2764;
  t2877 = t2862 + t2864 + t2874 + t2875;
  t3136 = -0.0265*t2604;
  t3137 = 0.2375*t2613;
  t3138 = t3136 + t3137;
  t2999 = 2.*t2998*t2972;
  t3000 = -1.*t2937*t2978;
  t3004 = -1.*t3003*t2981;
  t3007 = 0.0115*t3006;
  t3008 = -1.*t2925*t2937;
  t3010 = t3008 + t3009;
  t3011 = 0.0265*t3010;
  t3012 = t3000 + t3004 + t3007 + t3011;
  t3013 = 2.*t2956*t3012;
  t3026 = 2.*t2998*t2956;
  t3027 = 2.*t2987*t2972;
  t3028 = t3026 + t3027;
  t3046 = 2.*t3045*t2972;
  t3050 = 2.*t2956*t3049;
  t3051 = t3046 + t3050;
  t3056 = t744*t3040;
  t3165 = -0.0695*t2908;
  t3166 = 0.0265*t2912;
  t3167 = t3165 + t3166;
  t3055 = 2.*t3045*t2956;
  t3057 = t2522*t3043;
  t3066 = t3056 + t3057 + t3058 + t3059 + t3062 + t3064;
  t3067 = 2.*t2956*t3066;
  t3068 = 2.*t2972*t3049;
  t3070 = 2.*t3049*t3022;
  t3071 = t3055 + t3067 + t3068 + t3070;
  t3187 = -0.2375*t2925;
  t3188 = 0.0265*t2931;
  t3189 = t3187 + t3188;
  t3156 = 2.*t3045*t2987;
  t3157 = 2.*t2998*t3049;
  t3158 = t3156 + t2999 + t3013 + t3157;
  t2988 = 2.*t2956*t2987;
  t3023 = 2.*t2998*t3022;
  t3024 = t2988 + t2999 + t3013 + t3023;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=-0.25*Power(t2798,2)*t2803*var1[9] + 0.5*t2878*(2.*t2735*t2794 + 2.*Power(t2794,2) + 2.*t2801 + 2.*t2764*(t2618 + t2652 + t2670 + t2693 + t2522*t2711 - 1.*t2705*t744))*var1[9] - 0.25*t2766*t2798*t2803*var1[10] + 0.5*t2878*t2893*var1[10] - 0.25*t2798*t2803*t2858*var1[11] + 0.5*t2877*t2878*var1[11];
  p_output1[5]=-0.25*t3029*Power(t3034,2)*var1[9] + 0.5*t2975*(2.*t2957 + 2.*t2972*t3022 + 2.*Power(t3022,2) + 2.*t2956*(t2522*t2914 + t3058 + t3059 + t3062 + t3064 - 1.*t2920*t744))*var1[9] - 0.25*t3029*t3034*t3051*var1[12] + 0.5*t2975*t3071*var1[12] + 0.5*t2975*t3024*var1[13] - 0.25*t3028*t3029*t3034*var1[13];
  p_output1[6]=-0.25*t2766*t2798*t2803*var1[9] + 0.5*t2878*t2893*var1[9] - 0.25*Power(t2766,2)*t2803*var1[10] + 0.5*t2878*(2.*Power(t2695,2) + 2.*Power(t2743,2) + 2.*t2735*(t2770 + t2772 + t2783 + t2792 - 1.*t2522*t3103 - 1.*t2561*t744) + 2.*t2764*(t2585 + t2618 + t2652 + t2670 + t2693 + t3103*t744))*var1[10] - 0.25*t2766*t2803*t2858*var1[11] + 0.5*t2878*t3122*var1[11];
  p_output1[7]=0;
  p_output1[8]=-0.25*t2798*t2803*t2858*var1[9] + 0.5*t2877*t2878*var1[9] - 0.25*t2766*t2803*t2858*var1[10] + 0.5*t2878*t3122*var1[10] - 0.25*t2803*Power(t2858,2)*var1[11] + 0.5*t2878*(2.*Power(t2842,2) + 2.*Power(t2856,2) + 2.*t2735*(0.0265*t2836 + t3113 + 0.0115*(t2761 + t3115) - 1.*t2721*t3138) + 2.*t2764*(t2831 + 0.0115*(t2776 + t2832) + 0.0265*t2850 - 1.*t2754*t3138))*var1[11];
  p_output1[9]=0;
  p_output1[10]=0;
  p_output1[11]=-0.25*t3029*t3034*t3051*var1[9] + 0.5*t2975*t3071*var1[9] - 0.25*t3029*Power(t3051,2)*var1[12] + 0.5*t2975*(2.*Power(t3045,2) + 2.*Power(t3049,2) + 2.*t2956*(t3056 + t3058 + t3059 + t3062 + t3064 - 1.*t2522*t3167) + 2.*t2972*(t3016 + t3017 + t3018 + t3021 + t3047 - 1.*t3167*t744))*var1[12] - 0.25*t3028*t3029*t3051*var1[13] + 0.5*t2975*t3158*var1[13];
  p_output1[12]=0;
  p_output1[13]=0.5*t2975*t3024*var1[9] - 0.25*t3028*t3029*t3034*var1[9] - 0.25*t3028*t3029*t3051*var1[12] + 0.5*t2975*t3158*var1[12] - 0.25*Power(t3028,2)*t3029*var1[13] + 0.5*t2975*(2.*Power(t2987,2) + 2.*Power(t2998,2) + 2.*t2956*(0.0115*t2995 + t3000 + 0.0265*(t2969 + t3008) - 1.*t2924*t3189) + 2.*t2972*(0.0115*t2985 + t2989 + 0.0265*(t2991 + t3060) - 1.*t2962*t3189))*var1[13];
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
