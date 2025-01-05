  
% add path 
addpath('D:\Working\packages\calculation\frost-dev'); % change this path to where you install FROST
frost_addpath;
urdf = strcat('Robot_Assembly_v3_straight_leg.urdf'); % change this path to the path to urdf
export_path  = 'Expressions'; % change this path to where you want to store the expressions

robotModel = walker_class(urdf);

%%%% eom etc
robotModel.configureDynamics('DelayCoriolisSet', true);
robotModel.compile(export_path);

GCvec = SymExpression(zeros(7,1));
for i = 1:127
    GCvec = GCvec + robotModel.Fvec{ i};
end
GCvecFun = SymFunction('CGvec_five_link_walker',GCvec, robotModel.States.x, robotModel.States.dx); % change the name to what you want
GCvecFun.export(export_path);

expr = GenExpr5LinkWalker(robotModel,export_path);
expr.generate();


