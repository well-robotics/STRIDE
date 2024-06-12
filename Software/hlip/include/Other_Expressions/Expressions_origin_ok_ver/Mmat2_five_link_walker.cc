/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:36:39 GMT-05:00
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
  double t10;
  double t5;
  double t7;
  double t13;
  double t19;
  double t20;
  double t26;
  double t8;
  double t14;
  double t15;
  double t45;
  double t46;
  double t47;
  double t39;
  double t40;
  double t41;
  double t32;
  double t33;
  double t34;
  double t35;
  double t37;
  double t38;
  double t27;
  double t30;
  double t43;
  double t48;
  double t49;
  double t52;
  double t53;
  double t59;
  double t51;
  double t60;
  double t63;
  double t71;
  double t72;
  double t73;
  double t64;
  double t65;
  double t66;
  double t74;
  double t75;
  double t78;
  double t86;
  double t87;
  double t90;
  t10 = Cos(var1[2]);
  t5 = Cos(var1[3]);
  t7 = Sin(var1[2]);
  t13 = Sin(var1[3]);
  t19 = t10*t5;
  t20 = -1.*t7*t13;
  t26 = t19 + t20;
  t8 = t5*t7;
  t14 = t10*t13;
  t15 = t8 + t14;
  t45 = -0.00102*t5;
  t46 = -0.078722*t13;
  t47 = t45 + t46;
  t39 = -0.078722*t5;
  t40 = 0.00102*t13;
  t41 = t39 + t40;
  t32 = -1.*t5*t7;
  t33 = -1.*t10*t13;
  t34 = t32 + t33;
  t35 = 0.64654*t34*t26;
  t37 = 0.64654*t15*t26;
  t38 = t35 + t37;
  t27 = Power(t26,2);
  t30 = 0.64654*t27;
  t43 = t5*t41;
  t48 = t47*t13;
  t49 = t43 + t48;
  t52 = -1.*t5*t47;
  t53 = t41*t13;
  t59 = t52 + t53;
  t51 = 0.64654*t26*t49;
  t60 = 0.64654*t15*t59;
  t63 = t51 + t60;
  t71 = 0.64654*t34*t49;
  t72 = 0.64654*t26*t59;
  t73 = t71 + t72;
  t64 = 0.0006594708000000001*t15;
  t65 = -0.05089692188*t26;
  t66 = t64 + t65;
  t74 = -0.05089692188*t34;
  t75 = 0.0006594708000000001*t26;
  t78 = t74 + t75;
  t86 = -0.05089692188*t49;
  t87 = 0.0006594708000000001*t59;
  t90 = 0.000137 + t86 + t87;
  p_output1[0]=0.64654*Power(t15,2) + t30;
  p_output1[1]=t38;
  p_output1[2]=t63;
  p_output1[3]=t66;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t38;
  p_output1[8]=t30 + 0.64654*Power(t34,2);
  p_output1[9]=t73;
  p_output1[10]=t78;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t63;
  p_output1[15]=t73;
  p_output1[16]=0.000137 + 0.64654*Power(t49,2) + 0.64654*Power(t59,2);
  p_output1[17]=t90;
  p_output1[18]=0;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t66;
  p_output1[22]=t78;
  p_output1[23]=t90;
  p_output1[24]=0.00414438014445336;
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

#include "Mmat2_five_link_walker.hh"

namespace SymFunction
{

void Mmat2_five_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
