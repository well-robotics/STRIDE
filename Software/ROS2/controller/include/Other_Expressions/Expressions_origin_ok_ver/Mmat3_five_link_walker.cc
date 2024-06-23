/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:36:42 GMT-05:00
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
  double t8;
  double t14;
  double t18;
  double t31;
  double t34;
  double t35;
  double t15;
  double t19;
  double t26;
  double t60;
  double t63;
  double t64;
  double t49;
  double t51;
  double t53;
  double t40;
  double t41;
  double t43;
  double t45;
  double t46;
  double t48;
  double t37;
  double t38;
  double t59;
  double t65;
  double t66;
  double t68;
  double t70;
  double t71;
  double t67;
  double t72;
  double t73;
  double t83;
  double t84;
  double t86;
  double t74;
  double t75;
  double t78;
  double t87;
  double t90;
  double t91;
  double t99;
  double t100;
  double t101;
  t16 = Cos(var1[2]);
  t8 = Cos(var1[5]);
  t14 = Sin(var1[2]);
  t18 = Sin(var1[5]);
  t31 = t16*t8;
  t34 = -1.*t14*t18;
  t35 = t31 + t34;
  t15 = t8*t14;
  t19 = t16*t18;
  t26 = t15 + t19;
  t60 = -0.001112*t8;
  t63 = -0.078865*t18;
  t64 = t60 + t63;
  t49 = -0.078865*t8;
  t51 = 0.001112*t18;
  t53 = t49 + t51;
  t40 = -1.*t8*t14;
  t41 = -1.*t16*t18;
  t43 = t40 + t41;
  t45 = 0.64429*t43*t35;
  t46 = 0.64429*t26*t35;
  t48 = t45 + t46;
  t37 = Power(t35,2);
  t38 = 0.64429*t37;
  t59 = t8*t53;
  t65 = t64*t18;
  t66 = t59 + t65;
  t68 = -1.*t8*t64;
  t70 = t53*t18;
  t71 = t68 + t70;
  t67 = 0.64429*t35*t66;
  t72 = 0.64429*t26*t71;
  t73 = t67 + t72;
  t83 = 0.64429*t43*t66;
  t84 = 0.64429*t35*t71;
  t86 = t83 + t84;
  t74 = 0.00071645048*t26;
  t75 = -0.050811930850000006*t35;
  t78 = t74 + t75;
  t87 = -0.050811930850000006*t43;
  t90 = 0.00071645048*t35;
  t91 = t87 + t90;
  t99 = -0.050811930850000006*t66;
  t100 = 0.00071645048*t71;
  t101 = 0.000137 + t99 + t100;
  p_output1[0]=0.64429*Power(t26,2) + t38;
  p_output1[1]=t48;
  p_output1[2]=t73;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t78;
  p_output1[6]=0;
  p_output1[7]=t48;
  p_output1[8]=t38 + 0.64429*Power(t43,2);
  p_output1[9]=t86;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t91;
  p_output1[13]=0;
  p_output1[14]=t73;
  p_output1[15]=t86;
  p_output1[16]=0.000137 + 0.64429*Power(t66,2) + 0.64429*Power(t71,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t101;
  p_output1[20]=0;
  p_output1[21]=0;
  p_output1[22]=0;
  p_output1[23]=0;
  p_output1[24]=0;
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
  p_output1[35]=t78;
  p_output1[36]=t91;
  p_output1[37]=t101;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.00414507961941901;
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

#include "Mmat3_five_link_walker.hh"

namespace SymFunction
{

void Mmat3_five_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
