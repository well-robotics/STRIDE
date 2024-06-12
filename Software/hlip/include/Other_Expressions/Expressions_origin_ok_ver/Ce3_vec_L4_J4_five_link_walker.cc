/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:55 GMT-05:00
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
  double t1123;
  double t1136;
  double t1173;
  double t1166;
  double t1144;
  double t1145;
  double t1147;
  double t1152;
  double t1156;
  double t1157;
  double t1127;
  double t1137;
  double t1139;
  double t1182;
  double t1164;
  double t1184;
  double t1189;
  double t1193;
  double t1200;
  double t1202;
  double t1203;
  double t1204;
  double t1206;
  double t1208;
  double t1141;
  double t1159;
  double t1161;
  double t1167;
  double t1174;
  double t1179;
  double t1219;
  double t1209;
  double t1220;
  double t1240;
  double t1196;
  double t1258;
  double t1259;
  double t1261;
  double t1276;
  double t1279;
  double t1280;
  double t1283;
  double t1284;
  double t1285;
  double t1287;
  double t1288;
  double t1289;
  double t1313;
  double t1314;
  double t1315;
  double t1310;
  double t1311;
  double t1312;
  double t1316;
  double t1317;
  double t1319;
  double t1320;
  double t1321;
  double t1322;
  double t1323;
  double t1180;
  double t1197;
  double t1257;
  double t1262;
  double t1263;
  double t1264;
  double t1267;
  double t1268;
  double t1269;
  double t1270;
  double t1271;
  double t1272;
  double t1281;
  double t1291;
  double t1297;
  double t1354;
  double t1355;
  double t1356;
  double t1350;
  double t1351;
  double t1352;
  double t1301;
  t1123 = Cos(var1[4]);
  t1136 = Sin(var1[4]);
  t1173 = Cos(var1[3]);
  t1166 = Sin(var1[3]);
  t1144 = -1.*t1123;
  t1145 = 1. + t1144;
  t1147 = -0.16*t1145;
  t1152 = -0.167371*t1123;
  t1156 = 0.022663*t1136;
  t1157 = t1147 + t1152 + t1156;
  t1127 = -0.022663*t1123;
  t1137 = -0.007370999999999989*t1136;
  t1139 = t1127 + t1137;
  t1182 = Cos(var1[2]);
  t1164 = Sin(var1[2]);
  t1184 = t1173*t1123;
  t1189 = -1.*t1166*t1136;
  t1193 = t1184 + t1189;
  t1200 = t1123*t1157;
  t1202 = t1139*t1136;
  t1203 = t1200 + t1202;
  t1204 = -1.*t1123*t1166;
  t1206 = -1.*t1173*t1136;
  t1208 = t1204 + t1206;
  t1141 = -1.*t1123*t1139;
  t1159 = t1157*t1136;
  t1161 = t1141 + t1159;
  t1167 = t1123*t1166;
  t1174 = t1173*t1136;
  t1179 = t1167 + t1174;
  t1219 = -1.*t1164*t1193;
  t1209 = t1182*t1208;
  t1220 = t1209 + t1219;
  t1240 = -1.*t1164*t1208;
  t1196 = t1182*t1193;
  t1258 = -1.*t1173*t1123;
  t1259 = t1166*t1136;
  t1261 = t1258 + t1259;
  t1276 = -1.*t1166*t1139;
  t1279 = t1173*t1157;
  t1280 = t1276 + t1279;
  t1283 = -1.*t1173*t1139;
  t1284 = -1.*t1166*t1157;
  t1285 = t1283 + t1284;
  t1287 = t1173*t1139;
  t1288 = t1166*t1157;
  t1289 = t1287 + t1288;
  t1313 = 0.022663*t1123;
  t1314 = 0.007370999999999989*t1136;
  t1315 = t1313 + t1314;
  t1310 = -0.007370999999999989*t1123;
  t1311 = t1310 + t1156;
  t1312 = -1.*t1123*t1311;
  t1316 = t1315*t1136;
  t1317 = t1200 + t1312 + t1202 + t1316;
  t1319 = t1123*t1139;
  t1320 = t1123*t1315;
  t1321 = -1.*t1157*t1136;
  t1322 = t1311*t1136;
  t1323 = t1319 + t1320 + t1321 + t1322;
  t1180 = -1.*t1164*t1179;
  t1197 = t1180 + t1196;
  t1257 = 0.14994*t1161*t1220;
  t1262 = t1182*t1261;
  t1263 = t1240 + t1262;
  t1264 = 0.14994*t1203*t1263;
  t1267 = t1164*t1208;
  t1268 = t1267 + t1196;
  t1269 = 0.14994*t1161*t1268;
  t1270 = t1164*t1261;
  t1271 = t1209 + t1270;
  t1272 = 0.14994*t1203*t1271;
  t1281 = t1280*t1208;
  t1291 = t1289*t1193;
  t1297 = -1.*t1289*t1208;
  t1354 = t1173*t1315;
  t1355 = -1.*t1166*t1311;
  t1356 = t1354 + t1355;
  t1350 = t1166*t1315;
  t1351 = t1173*t1311;
  t1352 = t1350 + t1351;
  t1301 = -1.*t1280*t1261;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.14994*t1161*t1197 + 0.14994*t1203*t1220)*var2[0] + 0.5*(0.14994*t1161*(-1.*t1179*t1182 + t1219) + 0.14994*t1203*(-1.*t1182*t1193 + t1240))*var2[1])*var2[3];
  p_output1[3]=(0.5*(t1269 + t1272)*var2[0] + 0.5*(t1257 + t1264)*var2[1] + 0.5*(0.14994*t1203*(t1179*t1280 + t1281 + t1193*t1285 + t1291) + 0.14994*t1161*(-1.*t1193*t1280 - 1.*t1208*t1285 + t1297 + t1301))*var2[2])*var2[3];
  p_output1[4]=var2[3]*(0.5*(t1269 + t1272 + 0.14994*(t1179*t1182 + t1164*t1193)*t1317 + 0.14994*t1268*t1323)*var2[0] + 0.5*(t1257 + t1264 + 0.14994*t1197*t1317 + 0.14994*t1220*t1323)*var2[1] + 0.5*(0.14994*(-1.*t1208*t1280 - 1.*t1193*t1289)*t1317 + 0.14994*(t1193*t1280 + t1179*t1289)*t1323 + 0.14994*t1203*(t1281 + t1291 + t1179*t1352 + t1193*t1356) + 0.14994*t1161*(t1297 + t1301 - 1.*t1193*t1352 - 1.*t1208*t1356))*var2[2] + 0.5*(0.29988*t1161*t1317 + 0.29988*t1203*t1323)*var2[3] + 0.5*(0.0033980902199999994*t1317 - 0.0011052077399999983*t1323)*var2[4]);
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

#include "Ce3_vec_L4_J4_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L4_J4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
