/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:57 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t1156;
  double t1137;
  double t1147;
  double t1159;
  double t1184;
  double t1127;
  double t1196;
  double t1197;
  double t1198;
  double t1209;
  double t1220;
  double t1223;
  double t1152;
  double t1161;
  double t1167;
  double t1234;
  double t1224;
  double t1241;
  double t1249;
  double t1200;
  double t1268;
  double t1270;
  double t1271;
  double t1302;
  double t1305;
  double t1306;
  double t1307;
  double t1308;
  double t1309;
  double t1310;
  double t1313;
  double t1314;
  double t1315;
  double t1316;
  double t1317;
  double t1320;
  double t1322;
  double t1323;
  double t1325;
  double t1326;
  double t1327;
  double t1346;
  double t1347;
  double t1348;
  double t1342;
  double t1343;
  double t1267;
  double t1273;
  double t1275;
  double t1279;
  double t1281;
  double t1282;
  double t1286;
  double t1287;
  double t1288;
  double t1291;
  double t1292;
  double t1294;
  double t1296;
  double t1300;
  double t1318;
  double t1328;
  double t1332;
  double t1367;
  double t1368;
  double t1369;
  double t1363;
  double t1364;
  double t1365;
  double t1334;
  t1156 = Cos(var1[3]);
  t1137 = Cos(var1[4]);
  t1147 = Sin(var1[3]);
  t1159 = Sin(var1[4]);
  t1184 = Cos(var1[2]);
  t1127 = Sin(var1[2]);
  t1196 = t1156*t1137;
  t1197 = -1.*t1147*t1159;
  t1198 = t1196 + t1197;
  t1209 = -1.*t1137*t1147;
  t1220 = -1.*t1156*t1159;
  t1223 = t1209 + t1220;
  t1152 = t1137*t1147;
  t1161 = t1156*t1159;
  t1167 = t1152 + t1161;
  t1234 = -1.*t1127*t1198;
  t1224 = t1184*t1223;
  t1241 = t1224 + t1234;
  t1249 = -1.*t1127*t1223;
  t1200 = t1184*t1198;
  t1268 = -1.*t1156*t1137;
  t1270 = t1147*t1159;
  t1271 = t1268 + t1270;
  t1302 = -0.022663*t1137;
  t1305 = -0.007370999999999989*t1159;
  t1306 = t1302 + t1305;
  t1307 = -1.*t1147*t1306;
  t1308 = -1.*t1137;
  t1309 = 1. + t1308;
  t1310 = -0.16*t1309;
  t1313 = -0.167371*t1137;
  t1314 = 0.022663*t1159;
  t1315 = t1310 + t1313 + t1314;
  t1316 = t1156*t1315;
  t1317 = t1307 + t1316;
  t1320 = -1.*t1156*t1306;
  t1322 = -1.*t1147*t1315;
  t1323 = t1320 + t1322;
  t1325 = t1156*t1306;
  t1326 = t1147*t1315;
  t1327 = t1325 + t1326;
  t1346 = 0.022663*t1137;
  t1347 = 0.007370999999999989*t1159;
  t1348 = t1346 + t1347;
  t1342 = -0.007370999999999989*t1137;
  t1343 = t1342 + t1314;
  t1267 = 0.0033980902199999994*t1241;
  t1273 = t1184*t1271;
  t1275 = t1249 + t1273;
  t1279 = -0.0011052077399999983*t1275;
  t1281 = t1267 + t1279;
  t1282 = 0.5*var2[1]*t1281;
  t1286 = t1127*t1223;
  t1287 = t1286 + t1200;
  t1288 = 0.0033980902199999994*t1287;
  t1291 = t1127*t1271;
  t1292 = t1224 + t1291;
  t1294 = -0.0011052077399999983*t1292;
  t1296 = t1288 + t1294;
  t1300 = 0.5*var2[0]*t1296;
  t1318 = t1317*t1223;
  t1328 = t1327*t1198;
  t1332 = -1.*t1327*t1223;
  t1367 = t1156*t1348;
  t1368 = -1.*t1147*t1343;
  t1369 = t1367 + t1368;
  t1363 = t1147*t1348;
  t1364 = t1156*t1343;
  t1365 = t1363 + t1364;
  t1334 = -1.*t1317*t1271;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.0033980902199999994*(-1.*t1127*t1167 + t1200) - 0.0011052077399999983*t1241)*var2[0] + 0.5*(0.0033980902199999994*(-1.*t1167*t1184 + t1234) - 0.0011052077399999983*(-1.*t1184*t1198 + t1249))*var2[1])*var2[4];
  p_output1[3]=(t1282 + t1300 + 0.5*(-0.0011052077399999983*(t1167*t1317 + t1318 + t1198*t1323 + t1328) + 0.0033980902199999994*(-1.*t1198*t1317 - 1.*t1223*t1323 + t1332 + t1334))*var2[2])*var2[4];
  p_output1[4]=(t1282 + t1300 + 0.5*(-0.0011052077399999983*(t1318 + t1328 + t1167*t1365 + t1198*t1369) + 0.0033980902199999994*(t1332 + t1334 - 1.*t1198*t1365 - 1.*t1223*t1369))*var2[2] + 0.5*(-0.0011052077399999983*(t1137*t1306 - 1.*t1159*t1315 + t1159*t1343 + t1137*t1348) + 0.0033980902199999994*(t1159*t1306 + t1137*t1315 - 1.*t1137*t1343 + t1159*t1348))*var2[3])*var2[4];
  p_output1[5]=0;
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

#include "Ce3_vec_L4_J5_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L4_J5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
