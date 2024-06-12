/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:07:19 GMT-05:00
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
  double t16;
  double t13;
  double t14;
  double t17;
  double t21;
  double t7;
  double t22;
  double t24;
  double t25;
  double t15;
  double t18;
  double t19;
  double t20;
  double t26;
  double t27;
  double t32;
  double t35;
  double t36;
  double t37;
  double t40;
  double t41;
  double t58;
  double t59;
  double t54;
  double t55;
  double t60;
  double t61;
  double t62;
  double t63;
  double t65;
  double t66;
  double t67;
  double t68;
  double t72;
  double t73;
  double t74;
  double t75;
  double t76;
  double t56;
  double t57;
  double t64;
  double t69;
  double t70;
  double t46;
  double t47;
  double t48;
  double t49;
  double t50;
  double t51;
  double t52;
  double t53;
  double t71;
  double t77;
  double t78;
  double t80;
  double t81;
  double t82;
  double t91;
  double t92;
  double t93;
  double t94;
  double t95;
  double t85;
  double t86;
  double t87;
  double t88;
  double t89;
  double t79;
  double t83;
  double t84;
  double t106;
  double t107;
  double t108;
  double t90;
  double t96;
  double t97;
  double t109;
  double t110;
  double t111;
  double t120;
  double t121;
  double t122;
  double t98;
  double t99;
  double t100;
  double t112;
  double t113;
  double t114;
  double t123;
  double t124;
  double t125;
  double t131;
  double t132;
  double t133;
  t16 = Cos(var1[3]);
  t13 = Cos(var1[4]);
  t14 = Sin(var1[3]);
  t17 = Sin(var1[4]);
  t21 = Cos(var1[2]);
  t7 = Sin(var1[2]);
  t22 = t16*t13;
  t24 = -1.*t14*t17;
  t25 = t22 + t24;
  t15 = t13*t14;
  t18 = t16*t17;
  t19 = t15 + t18;
  t20 = t7*t19;
  t26 = t21*t25;
  t27 = t20 + t26;
  t32 = -1.*t13*t14;
  t35 = -1.*t16*t17;
  t36 = t32 + t35;
  t37 = t21*t36;
  t40 = t7*t25;
  t41 = t37 + t40;
  t58 = -1.*t13;
  t59 = 1. + t58;
  t54 = -1.*t16;
  t55 = 1. + t54;
  t60 = -0.2375*t59;
  t61 = -0.312342*t13;
  t62 = 0.0013409999999999984*t17;
  t63 = t60 + t61 + t62;
  t65 = -0.0265*t59;
  t66 = -0.025159*t13;
  t67 = 0.07484200000000002*t17;
  t68 = t65 + t66 + t67;
  t72 = -0.0265*t55;
  t73 = -0.0695*t14;
  t74 = -1.*t14*t63;
  t75 = t16*t68;
  t76 = t72 + t73 + t74 + t75;
  t56 = -0.0695*t55;
  t57 = 0.0265*t14;
  t64 = t16*t63;
  t69 = t14*t68;
  t70 = t56 + t57 + t64 + t69;
  t46 = t21*t19;
  t47 = -1.*t7*t25;
  t48 = t46 + t47;
  t49 = 0.19605*t27*t48;
  t50 = -1.*t7*t36;
  t51 = t50 + t26;
  t52 = 0.19605*t51*t41;
  t53 = t49 + t52;
  t71 = -1.*t70*t19;
  t77 = -1.*t76*t25;
  t78 = t71 + t77;
  t80 = t76*t36;
  t81 = t70*t25;
  t82 = t80 + t81;
  t91 = 0.0265*t13;
  t92 = t13*t68;
  t93 = 0.0695*t17;
  t94 = t63*t17;
  t95 = t91 + t92 + t93 + t94;
  t85 = -0.0695*t13;
  t86 = -1.*t13*t63;
  t87 = 0.0265*t17;
  t88 = t68*t17;
  t89 = t85 + t86 + t87 + t88;
  t79 = 0.19605*t41*t78;
  t83 = 0.19605*t27*t82;
  t84 = t79 + t83;
  t106 = 0.19605*t51*t78;
  t107 = 0.19605*t48*t82;
  t108 = t106 + t107;
  t90 = 0.19605*t89*t27;
  t96 = 0.19605*t95*t41;
  t97 = t90 + t96;
  t109 = 0.19605*t95*t51;
  t110 = 0.19605*t89*t48;
  t111 = t109 + t110;
  t120 = 0.19605*t95*t78;
  t121 = 0.19605*t89*t82;
  t122 = -0.00011 + t120 + t121;
  t98 = 0.014672774100000004*t27;
  t99 = 0.0002629030499999997*t41;
  t100 = t98 + t99;
  t112 = 0.0002629030499999997*t51;
  t113 = 0.014672774100000004*t48;
  t114 = t112 + t113;
  t123 = 0.0002629030499999997*t78;
  t124 = 0.014672774100000004*t82;
  t125 = -0.00011 + t123 + t124;
  t131 = 0.0002629030499999997*t95;
  t132 = 0.014672774100000004*t89;
  t133 = 0.00011 + t131 + t132;
  p_output1[0]=0.19605*Power(t27,2) + 0.19605*Power(t41,2);
  p_output1[1]=t53;
  p_output1[2]=t84;
  p_output1[3]=t97;
  p_output1[4]=t100;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t53;
  p_output1[8]=0.19605*Power(t48,2) + 0.19605*Power(t51,2);
  p_output1[9]=t108;
  p_output1[10]=t111;
  p_output1[11]=t114;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t84;
  p_output1[15]=t108;
  p_output1[16]=0.00011 + 0.19605*Power(t78,2) + 0.19605*Power(t82,2);
  p_output1[17]=t122;
  p_output1[18]=t125;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t97;
  p_output1[22]=t111;
  p_output1[23]=t122;
  p_output1[24]=0.00011 + 0.19605*Power(t89,2) + 0.19605*Power(t95,2);
  p_output1[25]=t133;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t100;
  p_output1[29]=t114;
  p_output1[30]=t125;
  p_output1[31]=t133;
  p_output1[32]=0.0012084923121822506;
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

#include "Mmat3_Robot_Assembly_v3_straight_leg_TPU_foot.hh"

namespace SymFunction
{

void Mmat3_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
