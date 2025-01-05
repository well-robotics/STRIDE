/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:44 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t1259;
  double t1231;
  double t1244;
  double t1276;
  double t1250;
  double t1290;
  double t1292;
  double t1298;
  double t1307;
  double t1316;
  double t1343;
  double t1344;
  double t1345;
  double t1386;
  double t1389;
  double t1416;
  double t1417;
  double t1418;
  double t1420;
  double t1399;
  double t1401;
  double t1411;
  double t1412;
  double t1296;
  double t1324;
  double t1332;
  double t1335;
  double t1340;
  double t1346;
  double t1347;
  double t1352;
  double t1353;
  double t1354;
  double t1356;
  double t1357;
  double t1358;
  double t1361;
  double t1366;
  double t1369;
  double t1370;
  double t1373;
  double t1382;
  double t1413;
  double t1423;
  double t1427;
  double t1428;
  double t1429;
  double t1433;
  double t1434;
  double t1435;
  double t1449;
  double t1450;
  double t1444;
  double t1445;
  double t1446;
  t1259 = Cos(var1[2]);
  t1231 = Cos(var1[5]);
  t1244 = Sin(var1[2]);
  t1276 = Sin(var1[5]);
  t1250 = -1.*t1231*t1244;
  t1290 = -1.*t1259*t1276;
  t1292 = t1250 + t1290;
  t1298 = t1259*t1231;
  t1307 = -1.*t1244*t1276;
  t1316 = t1298 + t1307;
  t1343 = t1231*t1244;
  t1344 = t1259*t1276;
  t1345 = t1343 + t1344;
  t1386 = -1.*t1231;
  t1389 = 1. + t1386;
  t1416 = -0.0265*t1389;
  t1417 = -0.025413*t1231;
  t1418 = -0.08282*t1276;
  t1420 = t1416 + t1417 + t1418;
  t1399 = -0.0695*t1389;
  t1401 = -0.15232*t1231;
  t1411 = -0.0010869999999999977*t1276;
  t1412 = t1399 + t1401 + t1411;
  t1296 = -0.0571880382*t1292;
  t1324 = -0.0007505843699999984*t1316;
  t1332 = t1296 + t1324;
  t1335 = 0.5*var2[5]*t1332;
  t1340 = 1.38102*t1292*t1316;
  t1346 = 1.38102*t1345*t1316;
  t1347 = t1340 + t1346;
  t1352 = 0.5*var2[0]*t1347;
  t1353 = Power(t1292,2);
  t1354 = 0.69051*t1353;
  t1356 = 0.69051*t1292*t1345;
  t1357 = Power(t1316,2);
  t1358 = 0.69051*t1357;
  t1361 = -1.*t1259*t1231;
  t1366 = t1244*t1276;
  t1369 = t1361 + t1366;
  t1370 = 0.69051*t1316*t1369;
  t1373 = t1354 + t1356 + t1358 + t1370;
  t1382 = 0.5*var2[1]*t1373;
  t1413 = t1231*t1412;
  t1423 = t1420*t1276;
  t1427 = t1413 + t1423;
  t1428 = 0.69051*t1292*t1427;
  t1429 = -1.*t1231*t1420;
  t1433 = t1412*t1276;
  t1434 = t1429 + t1433;
  t1435 = 0.69051*t1316*t1434;
  t1449 = -0.08282*t1231;
  t1450 = t1449 + t1411;
  t1444 = -0.0010869999999999977*t1231;
  t1445 = 0.08282*t1276;
  t1446 = t1444 + t1445;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(t1335 + t1352 + t1382 + 0.5*(t1428 + t1435)*var2[2]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(t1335 + t1352 + t1382 + 0.5*(t1428 + t1435 + 0.69051*t1345*(t1413 + t1423 + t1276*t1446 - 1.*t1231*t1450) + 0.69051*t1316*(-1.*t1276*t1412 + t1231*t1420 + t1231*t1446 + t1276*t1450))*var2[2]);
  p_output1[6]=0;
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

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
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
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce3_vec_L4_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L4_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
