/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:08 GMT-05:00
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
  double t18;
  double t15;
  double t16;
  double t19;
  double t26;
  double t8;
  double t27;
  double t31;
  double t32;
  double t17;
  double t20;
  double t22;
  double t24;
  double t33;
  double t37;
  double t43;
  double t47;
  double t48;
  double t49;
  double t53;
  double t59;
  double t80;
  double t81;
  double t75;
  double t76;
  double t83;
  double t84;
  double t86;
  double t87;
  double t91;
  double t92;
  double t93;
  double t95;
  double t100;
  double t101;
  double t102;
  double t105;
  double t106;
  double t78;
  double t79;
  double t90;
  double t96;
  double t98;
  double t66;
  double t67;
  double t68;
  double t70;
  double t71;
  double t72;
  double t73;
  double t74;
  double t99;
  double t107;
  double t110;
  double t112;
  double t115;
  double t116;
  double t129;
  double t130;
  double t131;
  double t132;
  double t133;
  double t120;
  double t121;
  double t124;
  double t125;
  double t126;
  double t111;
  double t117;
  double t119;
  double t151;
  double t152;
  double t153;
  double t128;
  double t136;
  double t138;
  double t154;
  double t155;
  double t156;
  double t169;
  double t170;
  double t172;
  double t141;
  double t142;
  double t143;
  double t158;
  double t159;
  double t160;
  double t173;
  double t174;
  double t175;
  double t184;
  double t185;
  double t186;
  t18 = Cos(var1[3]);
  t15 = Cos(var1[4]);
  t16 = Sin(var1[3]);
  t19 = Sin(var1[4]);
  t26 = Cos(var1[2]);
  t8 = Sin(var1[2]);
  t27 = t18*t15;
  t31 = -1.*t16*t19;
  t32 = t27 + t31;
  t17 = t15*t16;
  t20 = t18*t19;
  t22 = t17 + t20;
  t24 = t8*t22;
  t33 = t26*t32;
  t37 = t24 + t33;
  t43 = -1.*t15*t16;
  t47 = -1.*t18*t19;
  t48 = t43 + t47;
  t49 = t26*t48;
  t53 = t8*t32;
  t59 = t49 + t53;
  t80 = -1.*t15;
  t81 = 1. + t80;
  t75 = -1.*t18;
  t76 = 1. + t75;
  t83 = -0.2375*t81;
  t84 = -0.314514*t15;
  t86 = 0.0012709999999999978*t19;
  t87 = t83 + t84 + t86;
  t91 = -0.0265*t81;
  t92 = -0.025229*t15;
  t93 = 0.07701400000000003*t19;
  t95 = t91 + t92 + t93;
  t100 = -0.0265*t76;
  t101 = -0.0695*t16;
  t102 = -1.*t16*t87;
  t105 = t18*t95;
  t106 = t100 + t101 + t102 + t105;
  t78 = -0.0695*t76;
  t79 = 0.0265*t16;
  t90 = t18*t87;
  t96 = t16*t95;
  t98 = t78 + t79 + t90 + t96;
  t66 = t26*t22;
  t67 = -1.*t8*t32;
  t68 = t66 + t67;
  t70 = 0.19964*t37*t68;
  t71 = -1.*t8*t48;
  t72 = t71 + t33;
  t73 = 0.19964*t72*t59;
  t74 = t70 + t73;
  t99 = -1.*t98*t22;
  t107 = -1.*t106*t32;
  t110 = t99 + t107;
  t112 = t106*t48;
  t115 = t98*t32;
  t116 = t112 + t115;
  t129 = 0.0265*t15;
  t130 = t15*t95;
  t131 = 0.0695*t19;
  t132 = t87*t19;
  t133 = t129 + t130 + t131 + t132;
  t120 = -0.0695*t15;
  t121 = -1.*t15*t87;
  t124 = 0.0265*t19;
  t125 = t95*t19;
  t126 = t120 + t121 + t124 + t125;
  t111 = 0.19964*t59*t110;
  t117 = 0.19964*t37*t116;
  t119 = t111 + t117;
  t151 = 0.19964*t72*t110;
  t152 = 0.19964*t68*t116;
  t153 = t151 + t152;
  t128 = 0.19964*t126*t37;
  t136 = 0.19964*t133*t59;
  t138 = t128 + t136;
  t154 = 0.19964*t133*t72;
  t155 = 0.19964*t126*t68;
  t156 = t154 + t155;
  t169 = 0.19964*t133*t110;
  t170 = 0.19964*t126*t116;
  t172 = -0.000106 + t169 + t170;
  t141 = 0.015375074960000006*t37;
  t142 = 0.0002537424399999996*t59;
  t143 = t141 + t142;
  t158 = 0.0002537424399999996*t72;
  t159 = 0.015375074960000006*t68;
  t160 = t158 + t159;
  t173 = 0.0002537424399999996*t110;
  t174 = 0.015375074960000006*t116;
  t175 = -0.000106 + t173 + t174;
  t184 = 0.0002537424399999996*t133;
  t185 = 0.015375074960000006*t126;
  t186 = 0.000106 + t184 + t185;
  p_output1[0]=0.19964*Power(t37,2) + 0.19964*Power(t59,2);
  p_output1[1]=t74;
  p_output1[2]=t119;
  p_output1[3]=t138;
  p_output1[4]=t143;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t74;
  p_output1[8]=0.19964*Power(t68,2) + 0.19964*Power(t72,2);
  p_output1[9]=t153;
  p_output1[10]=t156;
  p_output1[11]=t160;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t119;
  p_output1[15]=t153;
  p_output1[16]=0.000106 + 0.19964*Power(t110,2) + 0.19964*Power(t116,2);
  p_output1[17]=t172;
  p_output1[18]=t175;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t138;
  p_output1[22]=t156;
  p_output1[23]=t172;
  p_output1[24]=0.000106 + 0.19964*Power(t126,2) + 0.19964*Power(t133,2);
  p_output1[25]=t186;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t143;
  p_output1[29]=t160;
  p_output1[30]=t175;
  p_output1[31]=t186;
  p_output1[32]=0.0012904185296106806;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=0;
  p_output1[36]=0;
  p_output1[37]=0;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0;
  p_output1[41]=0;
  p_output1[42]=0;
  p_output1[43]=0;
  p_output1[44]=0;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=0;
  p_output1[48]=0;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "Mmat3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Mmat3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
