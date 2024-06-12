/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:07:18 GMT-05:00
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
  double t8;
  double t5;
  double t6;
  double t9;
  double t7;
  double t10;
  double t13;
  double t16;
  double t17;
  double t18;
  double t28;
  double t29;
  double t35;
  double t36;
  double t37;
  double t38;
  double t30;
  double t31;
  double t32;
  double t33;
  double t22;
  double t23;
  double t24;
  double t25;
  double t26;
  double t27;
  double t19;
  double t20;
  double t34;
  double t39;
  double t40;
  double t42;
  double t43;
  double t44;
  double t41;
  double t45;
  double t46;
  double t53;
  double t54;
  double t55;
  double t47;
  double t48;
  double t49;
  double t56;
  double t57;
  double t58;
  double t64;
  double t65;
  double t66;
  t8 = Cos(var1[2]);
  t5 = Cos(var1[3]);
  t6 = Sin(var1[2]);
  t9 = Sin(var1[3]);
  t7 = t5*t6;
  t10 = -1.*t8*t9;
  t13 = t7 + t10;
  t16 = t8*t5;
  t17 = t6*t9;
  t18 = t16 + t17;
  t28 = -1.*t5;
  t29 = 1. + t28;
  t35 = -0.0695*t29;
  t36 = -0.15232*t5;
  t37 = 0.0011329999999999986*t9;
  t38 = t35 + t36 + t37;
  t30 = -0.0265*t29;
  t31 = -0.025367*t5;
  t32 = 0.08282*t9;
  t33 = t30 + t31 + t32;
  t22 = 0.69051*t13*t18;
  t23 = -1.*t5*t6;
  t24 = t8*t9;
  t25 = t23 + t24;
  t26 = 0.69051*t25*t18;
  t27 = t22 + t26;
  t19 = Power(t18,2);
  t20 = 0.69051*t19;
  t34 = -1.*t5*t33;
  t39 = -1.*t38*t9;
  t40 = t34 + t39;
  t42 = t5*t38;
  t43 = -1.*t33*t9;
  t44 = t42 + t43;
  t41 = 0.69051*t13*t40;
  t45 = 0.69051*t18*t44;
  t46 = t41 + t45;
  t53 = 0.69051*t18*t40;
  t54 = 0.69051*t25*t44;
  t55 = t53 + t54;
  t47 = 0.0007823478299999989*t13;
  t48 = 0.0571880382*t18;
  t49 = t47 + t48;
  t56 = 0.0571880382*t25;
  t57 = 0.0007823478299999989*t18;
  t58 = t56 + t57;
  t64 = 0.0007823478299999989*t40;
  t65 = 0.0571880382*t44;
  t66 = -0.00025 + t64 + t65;
  p_output1[0]=0.69051*Power(t13,2) + t20;
  p_output1[1]=t27;
  p_output1[2]=t46;
  p_output1[3]=t49;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t27;
  p_output1[8]=t20 + 0.69051*Power(t25,2);
  p_output1[9]=t55;
  p_output1[10]=t58;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t46;
  p_output1[15]=t55;
  p_output1[16]=0.00025 + 0.69051*Power(t40,2) + 0.69051*Power(t44,2);
  p_output1[17]=t66;
  p_output1[18]=0;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t49;
  p_output1[22]=t58;
  p_output1[23]=t66;
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

#include "Mmat2_Robot_Assembly_v3_straight_leg_TPU_foot.hh"

namespace SymFunction
{

void Mmat2_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
