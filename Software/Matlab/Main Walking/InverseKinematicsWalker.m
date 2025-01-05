function q0 = InverseKinematicsWalker(comZ, q0, torsoAngle, leftFootOffSet)
%,motorTorqueLimits, miu, half_foot_length)
% statics, require the spring stiffness, damping is not needed
%INITIALPOS  try to solve the IK by solving a set of nonlinear equations\
% use varargin to add external force to get more spring deformation

options = optimoptions('fsolve','MaxIterations',4000,'MaxFunctionEvaluations',20000,'Algorithm','Levenberg-Marquardt','Display','off');

comPosX_fun = @(q) COMPosition(q)*[1;0;0]; 
comPosZ_fun = @(q) COMPosition(q)*[0;0;1]; 
torsoAngle_fun = @(q) q(3); 

leftFootPosX_fun = @(q) [1,0] * pLeftToe(q); 
leftFootPosZ_fun = @(q) [0,1] * pLeftToe(q); 
rightFootPosX_fun = @(q) [1,0] * pRightToe(q); 
rightFootPosZ_fun = @(q) [0,1] * pRightToe(q); 

setBaseOri = @(q) q(3);

x0 = q0;
% stack all functions to constrain the robot such that single solution
%       2 + 4 + 12 + 22 + 6 + 4(6) + 2  = 22+ 12+ 10
stacked_fun = @(q) [% 100*setBaseOri(q);
                    100*comPosX_fun(q); 
                    100*(comPosZ_fun(q) - comZ); 
                    10*(torsoAngle_fun(q) - torsoAngle); 
                    leftFootPosX_fun(q) - leftFootOffSet(1); 
                    rightFootPosX_fun(q);
                    100*(leftFootPosZ_fun(q)- leftFootOffSet(2));
                    100*rightFootPosZ_fun(q);
                    ];

[X,fval,exitflag] = fsolve(stacked_fun,x0,options);
q0 = X;
end