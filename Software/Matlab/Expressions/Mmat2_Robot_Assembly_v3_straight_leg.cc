/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:05 GMT-05:00
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
  double t10;
  double t5;
  double t7;
  double t13;
  double t8;
  double t14;
  double t15;
  double t18;
  double t19;
  double t20;
  double t38;
  double t39;
  double t47;
  double t48;
  double t49;
  double t51;
  double t40;
  double t41;
  double t43;
  double t45;
  double t27;
  double t30;
  double t31;
  double t32;
  double t33;
  double t37;
  double t22;
  double t24;
  double t46;
  double t52;
  double t53;
  double t60;
  double t63;
  double t64;
  double t59;
  double t65;
  double t66;
  double t74;
  double t75;
  double t76;
  double t67;
  double t68;
  double t70;
  double t78;
  double t79;
  double t80;
  double t90;
  double t91;
  double t92;
  t10 = Cos(var1[2]);
  t5 = Cos(var1[3]);
  t7 = Sin(var1[2]);
  t13 = Sin(var1[3]);
  t8 = t5*t7;
  t14 = -1.*t10*t13;
  t15 = t8 + t14;
  t18 = t10*t5;
  t19 = t7*t13;
  t20 = t18 + t19;
  t38 = -1.*t5;
  t39 = 1. + t38;
  t47 = -0.0695*t39;
  t48 = -0.15232*t5;
  t49 = 0.0011329999999999986*t13;
  t51 = t47 + t48 + t49;
  t40 = -0.0265*t39;
  t41 = -0.025367*t5;
  t43 = 0.08282*t13;
  t45 = t40 + t41 + t43;
  t27 = 0.69051*t15*t20;
  t30 = -1.*t5*t7;
  t31 = t10*t13;
  t32 = t30 + t31;
  t33 = 0.69051*t32*t20;
  t37 = t27 + t33;
  t22 = Power(t20,2);
  t24 = 0.69051*t22;
  t46 = -1.*t5*t45;
  t52 = -1.*t51*t13;
  t53 = t46 + t52;
  t60 = t5*t51;
  t63 = -1.*t45*t13;
  t64 = t60 + t63;
  t59 = 0.69051*t15*t53;
  t65 = 0.69051*t20*t64;
  t66 = t59 + t65;
  t74 = 0.69051*t20*t53;
  t75 = 0.69051*t32*t64;
  t76 = t74 + t75;
  t67 = 0.0007823478299999989*t15;
  t68 = 0.0571880382*t20;
  t70 = t67 + t68;
  t78 = 0.0571880382*t32;
  t79 = 0.0007823478299999989*t20;
  t80 = t78 + t79;
  t90 = 0.0007823478299999989*t53;
  t91 = 0.0571880382*t64;
  t92 = -0.00025 + t90 + t91;
  p_output1[0]=0.69051*Power(t15,2) + t24;
  p_output1[1]=t37;
  p_output1[2]=t66;
  p_output1[3]=t70;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t37;
  p_output1[8]=t24 + 0.69051*Power(t32,2);
  p_output1[9]=t76;
  p_output1[10]=t80;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t66;
  p_output1[15]=t76;
  p_output1[16]=0.00025 + 0.69051*Power(t53,2) + 0.69051*Power(t64,2);
  p_output1[17]=t92;
  p_output1[18]=0;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t70;
  p_output1[22]=t80;
  p_output1[23]=t92;
  p_output1[24]=0.004987199723815391;
  p_output1[25]=0;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=0;
  p_output1[29]=0;
  p_output1[30]=0;
  p_output1[31]=0;
  p_output1[32]=0;
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

#include "Mmat2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Mmat2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
