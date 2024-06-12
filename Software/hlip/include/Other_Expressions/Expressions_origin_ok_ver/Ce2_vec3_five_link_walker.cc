/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:34 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t973;
  double t960;
  double t966;
  double t978;
  double t983;
  double t959;
  double t989;
  double t998;
  double t999;
  double t1030;
  double t1031;
  double t1032;
  double t1033;
  double t1034;
  double t972;
  double t979;
  double t981;
  double t982;
  double t1000;
  double t1008;
  double t1055;
  double t1052;
  double t1053;
  double t1056;
  double t1060;
  double t1061;
  double t1062;
  double t1070;
  double t1071;
  double t1072;
  double t1073;
  double t1074;
  double t1054;
  double t1057;
  double t1058;
  double t1059;
  double t1063;
  double t1064;
  double t1011;
  double t1027;
  double t1028;
  double t1092;
  double t1093;
  double t1094;
  double t1042;
  double t1038;
  double t1039;
  double t1040;
  double t1041;
  double t1043;
  double t1066;
  double t1067;
  double t1068;
  double t1107;
  double t1108;
  double t1109;
  double t1082;
  double t1078;
  double t1079;
  double t1080;
  double t1081;
  double t1083;
  double t1096;
  double t1097;
  double t1098;
  double t1100;
  double t1101;
  double t1103;
  double t1104;
  double t1105;
  double t1111;
  double t1112;
  double t1113;
  double t1115;
  double t1116;
  double t1118;
  double t1119;
  double t1120;
  double t1173;
  double t1174;
  double t1175;
  double t1177;
  double t1178;
  double t1179;
  double t1193;
  double t1194;
  double t1195;
  double t1197;
  double t1198;
  double t1199;
  double t1133;
  double t1134;
  double t1135;
  double t1129;
  double t1130;
  double t1131;
  double t1045;
  double t1046;
  double t1047;
  double t1048;
  double t1035;
  double t1036;
  double t1037;
  double t1141;
  double t1142;
  double t1153;
  double t1154;
  double t1155;
  double t1149;
  double t1150;
  double t1151;
  double t1085;
  double t1086;
  double t1087;
  double t1088;
  double t1075;
  double t1076;
  double t1077;
  double t1161;
  double t1162;
  double t1095;
  double t1110;
  double t1124;
  double t1125;
  double t1126;
  double t1127;
  double t1128;
  double t1132;
  double t1136;
  double t1137;
  double t1138;
  double t1139;
  double t1140;
  double t1143;
  double t1144;
  double t1145;
  double t1146;
  double t1147;
  double t1148;
  double t1152;
  double t1156;
  double t1157;
  double t1158;
  double t1159;
  double t1160;
  double t1163;
  double t1164;
  double t1167;
  double t1168;
  double t1169;
  double t1170;
  double t1171;
  double t1176;
  double t1180;
  double t1181;
  double t1183;
  double t1184;
  double t1185;
  double t1187;
  double t1188;
  double t1189;
  double t1190;
  double t1191;
  double t1196;
  double t1200;
  double t1201;
  double t1203;
  double t1204;
  double t1205;
  double t1232;
  double t1233;
  double t1234;
  double t1235;
  double t1236;
  double t1237;
  double t1238;
  double t1239;
  double t1166;
  double t1172;
  double t1182;
  double t1186;
  double t1192;
  double t1202;
  double t1206;
  double t1207;
  double t1029;
  double t1044;
  double t1049;
  double t1050;
  double t1212;
  double t1213;
  double t1214;
  double t1215;
  double t1069;
  double t1084;
  double t1089;
  double t1090;
  double t1218;
  double t1219;
  double t1220;
  double t1221;
  t973 = Cos(var1[3]);
  t960 = Cos(var1[4]);
  t966 = Sin(var1[3]);
  t978 = Sin(var1[4]);
  t983 = Sin(var1[2]);
  t959 = Cos(var1[2]);
  t989 = t973*t960;
  t998 = -1.*t966*t978;
  t999 = t989 + t998;
  t1030 = -1.*t960;
  t1031 = 1. + t1030;
  t1032 = 0.5*t1031;
  t1033 = 0.671885*t960;
  t1034 = t1032 + t1033;
  t972 = -1.*t960*t966;
  t979 = -1.*t973*t978;
  t981 = t972 + t979;
  t982 = t959*t981;
  t1000 = -1.*t983*t999;
  t1008 = t982 + t1000;
  t1055 = Cos(var1[5]);
  t1052 = Cos(var1[6]);
  t1053 = Sin(var1[5]);
  t1056 = Sin(var1[6]);
  t1060 = t1055*t1052;
  t1061 = -1.*t1053*t1056;
  t1062 = t1060 + t1061;
  t1070 = -1.*t1052;
  t1071 = 1. + t1070;
  t1072 = 0.5*t1071;
  t1073 = 0.671885*t1052;
  t1074 = t1072 + t1073;
  t1054 = -1.*t1052*t1053;
  t1057 = -1.*t1055*t1056;
  t1058 = t1054 + t1057;
  t1059 = t959*t1058;
  t1063 = -1.*t983*t1062;
  t1064 = t1059 + t1063;
  t1011 = -1.*t973*t983;
  t1027 = -1.*t959*t966;
  t1028 = t1011 + t1027;
  t1092 = t959*t973;
  t1093 = -1.*t983*t966;
  t1094 = t1092 + t1093;
  t1042 = t959*t999;
  t1038 = t960*t966;
  t1039 = t973*t978;
  t1040 = t1038 + t1039;
  t1041 = -1.*t983*t1040;
  t1043 = t1041 + t1042;
  t1066 = -1.*t1055*t983;
  t1067 = -1.*t959*t1053;
  t1068 = t1066 + t1067;
  t1107 = t959*t1055;
  t1108 = -1.*t983*t1053;
  t1109 = t1107 + t1108;
  t1082 = t959*t1062;
  t1078 = t1052*t1053;
  t1079 = t1055*t1056;
  t1080 = t1078 + t1079;
  t1081 = -1.*t983*t1080;
  t1083 = t1081 + t1082;
  t1096 = t973*t983;
  t1097 = t959*t966;
  t1098 = t1096 + t1097;
  t1100 = t983*t981;
  t1101 = t1100 + t1042;
  t1103 = t959*t1040;
  t1104 = t983*t999;
  t1105 = t1103 + t1104;
  t1111 = t1055*t983;
  t1112 = t959*t1053;
  t1113 = t1111 + t1112;
  t1115 = t983*t1058;
  t1116 = t1115 + t1082;
  t1118 = t959*t1080;
  t1119 = t983*t1062;
  t1120 = t1118 + t1119;
  t1173 = t1034*t966;
  t1174 = 0.171885*t973*t978;
  t1175 = t1173 + t1174;
  t1177 = t973*t1034;
  t1178 = -0.171885*t966*t978;
  t1179 = t1177 + t1178;
  t1193 = t1074*t1053;
  t1194 = 0.171885*t1055*t1056;
  t1195 = t1193 + t1194;
  t1197 = t1055*t1074;
  t1198 = -0.171885*t1053*t1056;
  t1199 = t1197 + t1198;
  t1133 = -1.*t983*t981;
  t1134 = -1.*t959*t999;
  t1135 = t1133 + t1134;
  t1129 = -1.*t959*t973;
  t1130 = t983*t966;
  t1131 = t1129 + t1130;
  t1045 = t1034*t960;
  t1046 = Power(t978,2);
  t1047 = 0.171885*t1046;
  t1048 = t1045 + t1047;
  t1035 = t1034*t978;
  t1036 = -0.171885*t960*t978;
  t1037 = t1035 + t1036;
  t1141 = -1.*t959*t1040;
  t1142 = t1141 + t1000;
  t1153 = -1.*t983*t1058;
  t1154 = -1.*t959*t1062;
  t1155 = t1153 + t1154;
  t1149 = -1.*t959*t1055;
  t1150 = t983*t1053;
  t1151 = t1149 + t1150;
  t1085 = t1074*t1052;
  t1086 = Power(t1056,2);
  t1087 = 0.171885*t1086;
  t1088 = t1085 + t1087;
  t1075 = t1074*t1056;
  t1076 = -0.171885*t1052*t1056;
  t1077 = t1075 + t1076;
  t1161 = -1.*t959*t1080;
  t1162 = t1161 + t1063;
  t1095 = 6.8522*t1028*t1094;
  t1110 = 6.8522*t1068*t1109;
  t1124 = Power(t1028,2);
  t1125 = 3.4261*t1124;
  t1126 = 3.4261*t1028*t1098;
  t1127 = Power(t1094,2);
  t1128 = 3.4261*t1127;
  t1132 = 3.4261*t1094*t1131;
  t1136 = 0.85216*t1135*t1101;
  t1137 = Power(t1043,2);
  t1138 = 0.85216*t1137;
  t1139 = Power(t1008,2);
  t1140 = 0.85216*t1139;
  t1143 = 0.85216*t1142*t1105;
  t1144 = Power(t1068,2);
  t1145 = 3.4261*t1144;
  t1146 = 3.4261*t1068*t1113;
  t1147 = Power(t1109,2);
  t1148 = 3.4261*t1147;
  t1152 = 3.4261*t1109*t1151;
  t1156 = 0.85216*t1155*t1116;
  t1157 = Power(t1083,2);
  t1158 = 0.85216*t1157;
  t1159 = Power(t1064,2);
  t1160 = 0.85216*t1159;
  t1163 = 0.85216*t1162*t1120;
  t1164 = t1125 + t1126 + t1128 + t1132 + t1136 + t1138 + t1140 + t1143 + t1145 + t1146 + t1148 + t1152 + t1156 + t1158 + t1160 + t1163;
  t1167 = Power(t973,2);
  t1168 = 0.1494*t1167;
  t1169 = Power(t966,2);
  t1170 = 0.1494*t1169;
  t1171 = t1168 + t1170;
  t1176 = -1.*t1175*t999;
  t1180 = -1.*t981*t1179;
  t1181 = t1176 + t1180;
  t1183 = t1175*t1040;
  t1184 = t999*t1179;
  t1185 = t1183 + t1184;
  t1187 = Power(t1055,2);
  t1188 = 0.1494*t1187;
  t1189 = Power(t1053,2);
  t1190 = 0.1494*t1189;
  t1191 = t1188 + t1190;
  t1196 = -1.*t1195*t1062;
  t1200 = -1.*t1058*t1199;
  t1201 = t1196 + t1200;
  t1203 = t1195*t1080;
  t1204 = t1062*t1199;
  t1205 = t1203 + t1204;
  t1232 = -3.70591*t959;
  t1233 = 3.4261*t1131*t1171;
  t1234 = 0.85216*t1142*t1181;
  t1235 = 0.85216*t1135*t1185;
  t1236 = 3.4261*t1151*t1191;
  t1237 = 0.85216*t1162*t1201;
  t1238 = 0.85216*t1155*t1205;
  t1239 = t1232 + t1233 + t1234 + t1235 + t1236 + t1237 + t1238;
  t1166 = -3.70591*t983;
  t1172 = 3.4261*t1028*t1171;
  t1182 = 0.85216*t1043*t1181;
  t1186 = 0.85216*t1008*t1185;
  t1192 = 3.4261*t1068*t1191;
  t1202 = 0.85216*t1083*t1201;
  t1206 = 0.85216*t1064*t1205;
  t1207 = t1166 + t1172 + t1182 + t1186 + t1192 + t1202 + t1206;
  t1029 = 0.51185934*t1028;
  t1044 = 0.85216*t1037*t1043;
  t1049 = 0.85216*t1048*t1008;
  t1050 = t1029 + t1044 + t1049;
  t1212 = 0.51185934*t1131;
  t1213 = 0.85216*t1048*t1135;
  t1214 = 0.85216*t1037*t1142;
  t1215 = t1212 + t1213 + t1214;
  t1069 = 0.51185934*t1068;
  t1084 = 0.85216*t1077*t1083;
  t1089 = 0.85216*t1088*t1064;
  t1090 = t1069 + t1084 + t1089;
  t1218 = 0.51185934*t1151;
  t1219 = 0.85216*t1088*t1155;
  t1220 = 0.85216*t1077*t1162;
  t1221 = t1218 + t1219 + t1220;
  p_output1[0]=var2[2]*(-0.5*(t1095 + 6.8522*t1094*t1098 + 1.70432*t1008*t1101 + 1.70432*t1043*t1105 + t1110 + 6.8522*t1109*t1113 + 1.70432*t1064*t1116 + 1.70432*t1083*t1120)*var2[0] - 0.5*t1164*var2[1] - 0.5*t1207*var2[2] - 0.5*t1050*var2[3] - 0.0732367608*t1008*var2[4] - 0.5*t1090*var2[5] - 0.0732367608*t1064*var2[6]);
  p_output1[1]=var2[2]*(-0.5*t1164*var2[0] - 0.5*(t1095 + t1110 + 6.8522*t1028*t1131 + 1.70432*t1008*t1135 + 1.70432*t1043*t1142 + 6.8522*t1068*t1151 + 1.70432*t1064*t1155 + 1.70432*t1083*t1162)*var2[1] - 0.5*t1239*var2[2] - 0.5*t1215*var2[3] - 0.0732367608*t1135*var2[4] - 0.5*t1221*var2[5] - 0.0732367608*t1155*var2[6]);
  p_output1[2]=(-0.5*t1207*var2[0] - 0.5*t1239*var2[1])*var2[2];
  p_output1[3]=(-0.5*t1050*var2[0] - 0.5*t1215*var2[1])*var2[2];
  p_output1[4]=(-0.0732367608*t1008*var2[0] - 0.0732367608*t1135*var2[1])*var2[2];
  p_output1[5]=(-0.5*t1090*var2[0] - 0.5*t1221*var2[1])*var2[2];
  p_output1[6]=(-0.0732367608*t1064*var2[0] - 0.0732367608*t1155*var2[1])*var2[2];
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

#include "Ce2_vec3_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
