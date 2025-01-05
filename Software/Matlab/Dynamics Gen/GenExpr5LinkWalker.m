classdef GenExpr5LinkWalker < handle
    %GENERATEEXPRESSIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        q
        dq
        x
        export_path
        robot
        mass
        Ndof = 7;
        g = 9.81;
    end
    
    properties (Access = private)
        motorIdx =  4:7;
        baseIdx = 1:3;
        % passing expressions
        COMvel
        COMoffset2center
        J_lmf
        J_rmf
    end
    
    methods
        
        function obj = GenExpr5LinkWalker(robot,export_path)
            obj.q = robot.States.x;
            obj.dq = robot.States.dx;
            xanddx = horzcat(obj.q',obj.dq');
            xanddx = xanddx';
            obj.x = SymVariable(xanddx);
            obj.robot = robot;
            obj.export_path = export_path;
            obj.mass = robot.mass;
        end
        
        function generate(obj)
            obj.PelvisExpressions();
            obj.ComExpressions();
            obj.FootExpressions();
            obj.HipExpressions();
            obj.legExpressions();
            obj.KneeExpressions();
            
            %  obj.genGradient();
            %  obj.genGradientEOM();
        end
        
        function PelvisExpressions(obj)
            % pelvis orientation
            expr = {};
            pelvis_ori = obj.q({'BaseRotY'});
            expr{end+1} = SymFunction('pelvis_ori',pelvis_ori,obj.q);
            J_pelvis = jacobian(pelvis_ori, obj.q);
            expr{end+1} = SymFunction('J_pelvis_ori',J_pelvis,obj.q);
            
            pelvis_rotV = J_pelvis*obj.dq;
            expr{end+1} = SymFunction('pelvis_rotvel',pelvis_rotV,obj.q,obj.dq);
            
            dJ_pelvis = jacobian(pelvis_rotV, obj.q);
            expr{end+1} = SymFunction('dJ_pelvis_ori',dJ_pelvis,obj.q,obj.dq);
            % pelvis position
            pelvis_pos = obj.q({'BasePosX','BasePosZ'});
            expr{end+1} = SymFunction('pelvis_pos',pelvis_pos,obj.q);
            
            J_pelvis_pos = jacobian(pelvis_pos, obj.q);
            expr{end+1} = SymFunction('J_pelvis_pos',J_pelvis_pos,obj.q);
            
            pelvis_vel = J_pelvis_pos*obj.dq;
            expr{end+1} = SymFunction('pelvis_vel',pelvis_vel,obj.q,obj.dq);
            
            dJ_pelvis_pos = jacobian(pelvis_vel,obj.q);
            expr{end+1} = SymFunction('dJ_pelvis_pos',dJ_pelvis_pos,obj.q,obj.dq);
            
            obj.processExport(expr, obj.export_path);
        end
        
        function ComExpressions(obj)
            expr = {};
            
            r = obj.robot; X = obj.q; dX = obj.dq;
            % compute com position and jacobian in base frame
            COM_pos = r.getComPosition; %1x3
            J_COM = jacobian(COM_pos,X);
            
            COM_vel = J_COM*dX;
            obj.COMvel = COM_vel;
            dJ_COM = jacobian(COM_vel,X);
            
            expr{end+1} = SymFunction('COMPosition',COM_pos,X);
            expr{end+1}  = SymFunction('J_COMPosition',J_COM,X);
            
            expr{end+1}  = SymFunction('COM_velocity',COM_vel,X,dX);
            expr{end+1}  = SymFunction('dJ_COMPosition',dJ_COM,X,dX);
            
            obj.processExport(expr, obj.export_path);
        end
        
        function FootExpressions(obj)
            expr = {};
            r = obj.robot;
            
            %% foot heel and toes
            Left_toe = r.getCartesianPosition(r.ContactPoints.LeftToe);
            Right_toe = r.getCartesianPosition(r.ContactPoints.RightToe);
            expr{end+1} = SymFunction('pLeftToe', Left_toe(1:2:3), obj.q);
            expr{end+1} = SymFunction('pRightToe', Right_toe(1:2:3), obj.q);
            
            J_Left_toe = jacobian(Left_toe(1:2:3), obj.q);
            J_Right_toe = jacobian(Right_toe(1:2:3), obj.q);
            expr{end+1} = SymFunction('J_leftToe', J_Left_toe, obj.q);
            expr{end+1} = SymFunction('J_rightToe', J_Right_toe, obj.q);
            
            v_Left_toe = J_Left_toe*obj.dq;
            v_Right_toe = J_Right_toe*obj.dq;
            expr{end+1} = SymFunction('vLeftToe',v_Left_toe,obj.x);
            expr{end+1} = SymFunction('vRightToe',v_Right_toe,obj.x);
            
            dJ_Left_toe = jacobian(v_Left_toe , obj.q);
            dJ_Right_toe = jacobian(v_Right_toe, obj.q);
            expr{end+1} = SymFunction('dJ_leftToe',dJ_Left_toe,obj.x);
            expr{end+1} = SymFunction('dJ_rightToe',dJ_Right_toe,obj.x);
            
            obj.processExport(expr, obj.export_path);
            
        end
        
        function KneeExpressions(obj)
            expr = {};
            r = obj.robot;
            
            %% foot heel and toes
            Left_knee = r.getCartesianPosition(r.ContactPoints.LeftKnee);
            Right_knee = r.getCartesianPosition(r.ContactPoints.RightKnee);
            expr{end+1} = SymFunction('pLeftKnee', Left_knee, obj.q);
            expr{end+1} = SymFunction('pRightKnee', Right_knee, obj.q);
            
            TorsoTop = r.getCartesianPosition(r.ContactPoints.Torso);
            expr{end+1} = SymFunction('pTorsoTop', TorsoTop, obj.q);
            
            obj.processExport(expr, obj.export_path);
        end
        
        function HipExpressions(obj)
            expr = {};
            RobotModel = obj.robot;
            
            hipPos = RobotModel.getCartesianPosition(RobotModel.ContactPoints.HipPos);
            expr{end+1} = SymFunction('pHip', hipPos, obj.q);
            
            J_hipPos = jacobian(hipPos, obj.q);
            expr{end+1} = SymFunction('J_hip', J_hipPos, obj.q);
            
            v_hipPos = J_hipPos*obj.dq;
            expr{end+1} = SymFunction('vHip', v_hipPos, obj.x);
            
            dJ_hipPos = jacobian(v_hipPos, obj.q);
            expr{end+1} = SymFunction('dJ_hip', dJ_hipPos, obj.x);
            
            obj.processExport(expr, obj.export_path);
            
        end
        
        function legExpressions(obj)
            expr = {};
            
            RobotModel = obj.robot;
            Left_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.LeftToe);
            Right_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.RightToe);
            hipPos = RobotModel.getCartesianPosition(RobotModel.ContactPoints.HipPos);
            
            LegL = SymExpression(zeros(2,1));
            temp1 = hipPos - Right_toe;
            temp2 = hipPos - Left_toe;
            LegL(1,1) = sqrt(temp1(1,1).^2 + temp1(1,2).^2 + temp1(1,3).^2);
            LegL(2,1) = sqrt(temp2(1,1).^2 + temp2(1,2).^2 + temp2(1,3).^2);
            
            expr{end+1} = SymFunction('LegL', LegL, obj.q);
            
            J_LegL = jacobian(LegL, obj.q);
            expr{end+1} = SymFunction('J_LegL', J_LegL, obj.q);
            
            v_LegL = J_LegL*obj.dq;
            expr{end+1} = SymFunction('vLegL', v_LegL, obj.x);
            
            dJ_LegL = jacobian(v_LegL, obj.q);
            expr{end+1} = SymFunction('dJ_LegL', dJ_LegL, obj.x);
            
            obj.processExport(expr, obj.export_path);
            
        end
        
        function genGradient(obj)
            %%% A(q) = inv(M)*J'*inv(J*inv(M)*J')*J*dq;
            %%% generate gradient of Aq;
            RobotModel = obj.robot;
            X = obj.q; dX = obj.dq;
            
            Left_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.LeftToe);
            J = jacobian(Left_toe,obj.q);
            v_Left_toe = J*obj.dq;
            dJ_Left_toe = jacobian(v_Left_toe, obj.q);
            
            JJ = J(1:2:3,:); %%%%select the x and z components \
            dJJ = dJ_Left_toe(1:2:3,:);
            M = RobotModel.Mmat;
            
            gradM = [];
            for i = 1:obj.Ndof
                gradM = [gradM; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end
            
            for i = 1:obj.Ndof
                for j = 1:obj.Ndof
                    gradM(j + (i-1)*obj.Ndof,: ) = jacobian(M(i,j), X);
                end
            end  %             gradM = gradientt(M, X);
            gradMfun = SymFunction('grad_Mq',gradM,obj.x);
            gradMfun.export(obj.export_path);
            %
            gradJ = [];
            for i = 1:2
                gradJ = [gradJ; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end    % gradJ = gradientt(JJ, X);
            for i = 1:2
                for j = 1:obj.Ndof
                    gradJ(j + (i-1)*obj.Ndof,: ) = jacobian(JJ(i,j), X);
                end
            end
            gradJfun = SymFunction('grad_J_left_q',gradJ,obj.x);
            gradJfun.export(obj.export_path);
            
            graddJ = [];
            for i = 1:2
                graddJ = [graddJ; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end    % gradJ = gradientt(JJ, X);
            for i = 1:2
                for j = 1:obj.Ndof
                    graddJ(j + (i-1)*obj.Ndof,: ) = jacobian(dJJ(i,j), X);
                end
            end
            graddJfun = SymFunction('grad_dJ_left_q',gradJ,obj.x);
            graddJfun.export(obj.export_path);
            
            graddJ_dq = [];
            for i = 1:2
                graddJ_dq = [graddJ_dq; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end    % gradJ = gradientt(JJ, X);
            for i = 1:2
                for j = 1:obj.Ndof
                    graddJ_dq(j + (i-1)*obj.Ndof,: ) = jacobian(dJJ(i,j), dX);
                end
            end
            graddJfun = SymFunction('grad_dJ_left_dq',graddJ_dq,obj.x);
            graddJfun.export(obj.export_path);
            
            %             invM = inv(M);
            %             invJMJ = inv(JJ*invM*JJ');
            %             Aq = invM*JJ'*invJMJ*JJ*dX;
            
            %             grad = jacobian(Aq, X);
            %             gradfun = SymFunction('grad_Aq',grad,obj.x);
            %             gradfun.export(obj.export_path);
            
            
            Right_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.RightToe);
            J = jacobian(Right_toe,obj.q);
            v_Right_toe = J*obj.dq;
            dJ_Right_toe = jacobian(v_Right_toe, obj.q);
            
            JJ = J(1:2:3,:); %%%%select the x and z components \
            dJJ = dJ_Right_toe(1:2:3,:);
            
            gradJ = [];
            for i = 1:2
                gradJ = [gradJ; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end    % gradJ = gradientt(JJ, X);
            for i = 1:2
                for j = 1:obj.Ndof
                    gradJ(j + (i-1)*obj.Ndof,: ) = jacobian(JJ(i,j), X);
                end
            end
            gradJfun = SymFunction('grad_J_right_q',gradJ,obj.x);
            gradJfun.export(obj.export_path);
            
            graddJ = [];
            for i = 1:2
                graddJ = [graddJ; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end    % gradJ = gradientt(JJ, X);
            for i = 1:2
                for j = 1:obj.Ndof
                    graddJ(j + (i-1)*obj.Ndof,: ) = jacobian(dJJ(i,j), X);
                end
            end
            graddJfun = SymFunction('grad_dJ_right_q',gradJ,obj.x);
            graddJfun.export(obj.export_path);
            
            graddJ_dq = [];
            for i = 1:2
                graddJ_dq = [graddJ_dq; SymExpression(zeros(obj.Ndof, obj.Ndof))];    % stack joint jacbians
            end    % gradJ = gradientt(JJ, X);
            for i = 1:2
                for j = 1:obj.Ndof
                    graddJ_dq(j + (i-1)*obj.Ndof,: ) = jacobian(dJJ(i,j), dX);
                end
            end
            graddJfun = SymFunction('grad_dJ_right_dq',graddJ_dq,obj.x);
            graddJfun.export(obj.export_path);
        end
        
        function genGradientEOM(obj)
            RobotModel = obj.robot;
            X = obj.q; dX = obj.dq;
            
            GCvec = SymExpression(zeros(obj.Ndof,1 ));
            for i = 1:22
                GCvec = GCvec - RobotModel.Fvec{i};
            end
            
            GCvec_Fun = SymFunction('CGvec_robot_assembly_v3', GCvec, obj.x);
            GCvec_Fun.export(obj.export_path);
            
            grad_CG_q = jacobian(GCvec, obj.q);
            grad_CG_qFun = SymFunction('grad_CG_q', grad_CG_q, obj.x);
            grad_CG_qFun.export(obj.export_path);
            grad_CG_dq = jacobian(GCvec, obj.dq);
            grad_CG_dqFun = SymFunction('grad_CG_dq', grad_CG_dq, obj.x);
            grad_CG_dqFun.export(obj.export_path);
            
        end
    end
    
    methods %%%% Discrete mechanics based terms
        
        function DiscreteMechanics(obj)
            expr = {};
            h = SymVariable('h'); %%%% discretization
            qk = SymVariable('qk', obj.Ndof);
            qkp1 = SymVariable('qkp1', obj.Ndof);
            
            %%% kinetic energy
            K = 1/2*transpose(obj.dq)*obj.robot.Mmat*obj.dq;
            %%% potential energy
            comPos = obj.robot.getComPosition;
            P = comPos(3)*obj.g*obj.mass;
            
            %%% Lagragain
            L = K - P;
            
            %%% discrete Lagraingain
            D2L = transpose(jacobian(L, obj.dq));
            
            Ld = h*L;
            Ld = Ld.subs(obj.q, (qk + qkp1)/2);
            Ld = Ld.subs(obj.dq, (qkp1 - qk)/h);
            
            D1Ld = transpose(jacobian(Ld, qk)); %%% vector
            D2Ld = transpose(jacobian(Ld, qkp1));
            
            d1_D1Ld = jacobian(D1Ld, qk);
            d2_D1Ld = jacobian(D1Ld, qkp1);
            
            d1_D2Ld = jacobian(D2Ld, qk);
            d2_D2Ld = jacobian(D2Ld, qkp1);
            
            expr{end+1} = SymFunction('D2L', D2L, {obj.q, obj.dq});
            
            expr{end+1} = SymFunction('D1Ld', D1Ld, {qk, qkp1, h});
            expr{end+1} = SymFunction('D2Ld', D2Ld, {qk, qkp1, h});
            
            expr{end+1} = SymFunction('d1_D1Ld', d1_D1Ld, {qk, qkp1, h});
            expr{end+1} = SymFunction('d2_D1Ld', d2_D1Ld, {qk, qkp1, h});
            
            expr{end+1} = SymFunction('d1_D2Ld', d1_D2Ld, {qk, qkp1, h});
            expr{end+1} = SymFunction('d2_D2Ld', d2_D2Ld, {qk, qkp1, h});
             %%% hessian matrix for D1Ld 
            for i = 1:length(qk) 
                Hi = jacobian(jacobian(D1Ld(i), qkp1), qkp1); 
                name = ['H', int2str(i), 'd2_D1Ld'];
                expr{end+1} = SymFunction(name, Hi, {qk, qkp1, h});
            end 
            
        
            path = [obj.export_path, '\+DM'];
            obj.processExport(expr, path);
            
%             addpath('C:\Users\Xiaobin\Dropbox\XiaobinNew\cassie_walking_stack\modules\cassie_description\MATLAB\symbolic');
%             %%% path in X\Onedrive is problematic 
%             path = 'C:\Users\Xiaobin\Desktop\temp';
%             outpath = [path, '\+Mfile']; 
%             for i = 1:length(expr)
%                  convertFile(expr{i}.Name, path, outpath)
%             end
        end
        
        function convertMex2Mfile(obj) 
            %%% convert some mex fill's original cc files to .m files 
            names = {};
           % names{end + 1} = 'D2Ld';
            names{end + 1} = 'J_leftToe';
            names{end + 1} = 'J_rightToe';
            
           
            path = 'C:\Users\Xiaobin\Desktop\temp';
            outpath = [path, '\+Mfile']; 
            for i = 1:length(names)
                 convertFile(names{i}, path, outpath)
            end
        end 
        
        function outputVec(obj)
            expr = {};
            r = obj.robot;
            
            COM_pos = r.getComPosition; %1x3
            pelvis_ori = obj.q({'BaseRotY'});
            Left_toe = r.getCartesianPosition(r.ContactPoints.LeftToe);
            Right_toe = r.getCartesianPosition(r.ContactPoints.RightToe);
            
            yLeftStance = [  Right_toe(1);
                Right_toe(3);
                COM_pos(3);
                pelvis_ori;];
            
            yRightStance = [  Left_toe(1);
                Left_toe(3);
                COM_pos(3);
                pelvis_ori;];
            
            JyLeft = jacobian(yLeftStance, obj.q); 
            JyRight = jacobian(yRightStance, obj.q); 
            
            expr{end+1} = SymFunction('yLeftStance', yLeftStance, obj.q);
            expr{end+1} = SymFunction('yRightStance', yRightStance, obj.q);
            
             expr{end+1} = SymFunction('JyLeftStance', JyLeft, obj.q);
            expr{end+1} = SymFunction('JyRightStance', JyRight, obj.q);
            
            
            path = [obj.export_path, '/+DM'];
            obj.processExport(expr, path);
            
        end
        
        function toeExpr(obj) 
            expr = {};
            r = obj.robot;
            
            %% foot heel and toes
            Left_toe = r.getCartesianPosition(r.ContactPoints.LeftToe);
            Right_toe = r.getCartesianPosition(r.ContactPoints.RightToe);
            
            expr{end+1} = SymFunction('pLeftToe', Left_toe(1:2:3), obj.q);
            expr{end+1} = SymFunction('pRightToe', Right_toe(1:2:3), obj.q);
            
            obj.processExport(expr, obj.export_path);
        end 
        
        function getHessian(obj)
            
            expr = {};
            r = obj.robot;
             COM_pos = r.getComPosition; %1x3

            %% foot heel and toes
            Left_toe = r.getCartesianPosition(r.ContactPoints.LeftToe);
            Right_toe = r.getCartesianPosition(r.ContactPoints.RightToe);
            Hi = jacobian(jacobian(Left_toe(1), obj.q), obj.q);
            expr{end+1} = SymFunction('H_leftToeX', Hi, obj.q);
            Hi = jacobian(jacobian(Left_toe(3), obj.q), obj.q);
            expr{end+1} = SymFunction('H_leftToeZ', Hi, obj.q);
            
            Hi = jacobian(jacobian(Right_toe(1), obj.q), obj.q);
            expr{end+1} = SymFunction('H_rightToeX', Hi, obj.q);
            Hi = jacobian(jacobian(Right_toe(3), obj.q), obj.q);
            expr{end+1} = SymFunction('H_rightToeZ', Hi, obj.q);
            
            Hi = jacobian(jacobian(COM_pos(3), obj.q), obj.q);
            expr{end+1} = SymFunction('H_COMZ', Hi, obj.q);
            pelvis_ori = obj.q({'BaseRotY'});

            obj.processExport(expr, [obj.export_path, '\+Hessian']); 
        end 
    end
    
    
    methods %%%%%% generation equation of motion in the normal coordinates
        
        function genNormalEOM(obj)
            RobotModel = obj.robot;
            X = obj.q;
            dX = obj.dq;
            Yvar = SymVariable('yvar', [obj.Ndof, 1]);
            dYvar = SymVariable('dyvar', [obj.Ndof,1]);
            
            %%%% normal coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% Y = [ zcom, leftFootXZ, rightFootXZ, qtorso, xcom];%%%%%
            COM_pos = obj.robot.getComPosition; %1x3
            pelvis_ori = obj.q({'BaseRotY'});
            Left_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.LeftToe);
            Right_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.RightToe);
            
            %             Y = [COM_pos(3);
            %                  Left_toe(1);
            %                  Left_toe(3);
            %                  Right_toe(1);
            %                  Right_toe(3);
            %                  pelvis_ori;
            %                  COM_pos(1)];
            
            Y = [COM_pos(1);
                COM_pos(3);
                pelvis_ori;
                Right_toe(1);
                Right_toe(3);
                Left_toe(1);
                Left_toe(3);
                ];
            
            JY = jacobian(Y, obj.q);
            
            dJY = SymExpression(zeros(obj.Ndof));
            for i = 1:obj.Ndof
                for j = 1: obj.Ndof
                    dJY(i,j) = jacobian(JY(i,j),obj.q)*obj.dq;
                end
            end
            
            % dJY = jacobian(JY*obj.dq, obj.q);
            dY = JY*obj.dq;
            A = inv(JY); %%%% dq = A*dY
            
            Adot = SymExpression(zeros(obj.Ndof));
            for i = 1:obj.Ndof
                for j = 1: obj.Ndof
                    Adot(i,j) = jacobian(A(i,j),obj.q)*obj.dq;
                end
            end
            
            % Adot = jacobian(A*obj.dq, obj.q); %%%
            
            Gnormal = SymExpression(zeros(obj.Ndof,1));
            Gnormal(2) = obj.robot.mass*obj.g;
            
            J_leftFoot = SymExpression(zeros(2, obj.Ndof));
            J_rightFoot = SymExpression(zeros(2, obj.Ndof));
            J_rightFoot(1, 4) = 1;
            J_rightFoot(2, 5) = 1;
            J_leftFoot(1, 6) = 1;
            J_leftFoot(2, 7) = 1;
            
            Mnormal = A'*obj.robot.Mmat*A;
            Bnormal = A'*obj.robot.Gmap.Control.u;
            
            %%%%%%%%%% partial M/ partial y = par M/par q * par q/par y
            %            gradM = SymExpression(zeros(obj.Ndof*obj.Ndof,obj.Ndof));
            %             for i = 2:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     gradM(j + (i-1)*obj.Ndof,: ) = jacobian(Mnormal(i,j), obj.q)*A;
            %                 end
            %             end  %             gradM = gradientt(M, X);
            %             TijkNormal = SymExpression(zeros( obj.Ndof*obj.Ndof, obj.Ndof));
            %             for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     for k = 1:obj.Ndof
            %                         TijkNormal(j + (i-1)*obj.Ndof, k) = 1/2*( gradM(j + (i-1)*obj.Ndof, k) + gradM(k + (i-1)*obj.Ndof, j) - gradM(k + (j-1)*obj.Ndof, i));
            %                     end
            %                 end
            %             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             corilisMat1 = SymExpression(zeros(obj.Ndof, obj.Ndof));
            %             corilisMat2 = SymExpression(zeros(obj.Ndof, obj.Ndof));
            %             corilisMat3 = SymExpression(zeros(obj.Ndof, obj.Ndof));
            %
            %             for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     corilisMat1(i,j) = 1/2* gradM(j + (i-1)*obj.Ndof, :)*dYvar;
            %                 end
            %             end
            %             for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     corilisMat2(i,j) = 1/2* gradM( (1:obj.Ndof) + (i-1)*obj.Ndof, j)*dYvar;
            %                 end
            %             end
            %              for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     corilisMat3(i,j) = -1/2*gradM((1:obj.Ndof)  + (j-1)*obj.Ndof, i)*dYvar;
            %                 end
            %             end
            %             corilisVec = SymExpression(zeros(obj.Ndof, 1));
            %             corilisVec  = corilisMat1*dYvar +  corilisMat2*dYvar +  corilisMat3*dYvar;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             CeSum = obj.CalculateNormalCoriolisRobotLink(Mnormal, A, dYvar);
            %             Ce1 = cell(obj.Ndof,1);
            %             Ce2 = cell(obj.Ndof,1);
            %             Ce3 = cell(obj.Ndof,1);
            %             Ce1sum = SymExpression(zeros(obj.Ndof, 1));
            %             Ce2sum = SymExpression(zeros(obj.Ndof, 1));
            %             Ce3sum = SymExpression(zeros(obj.Ndof, 1));
            %
            %             for k = 1:obj.Ndof
            %                 corilisVec1 = SymExpression(zeros(obj.Ndof, 1));
            %                 for i = 1:obj.Ndof
            %                     for j = 1:obj.Ndof
            %                         corilisVec1(i) = corilisVec1(i) + dYvar(j)*1/2* gradM(j + (i-1)*obj.Ndof, k)*dYvar(k);
            %                     end
            %                 end
            %                 Ce1{k} = corilisVec1;
            %             end
            %             for j = 1:obj.Ndof
            %                 Ce1sum = Ce1sum + Ce1{k};
            %             end
            %              for k = 1:obj.Ndof
            %                 corilisVec2 = SymExpression(zeros(obj.Ndof, 1));
            %                 for i = 1:obj.Ndof
            %                     for j = 1:obj.Ndof
            %                         corilisVec2(i) = corilisVec2(i) + dYvar(j)*1/2* gradM(k + (i-1)*obj.Ndof, j)*dYvar(k);
            %                     end
            %                 end
            %                 Ce2{k} = corilisVec2;
            %              end
            %             for j = 1:obj.Ndof
            %                 Ce2sum = Ce2sum + Ce2{k};
            %             end
            %
            %             for k = 1:obj.Ndof
            %                 corilisVec3 = SymExpression(zeros(obj.Ndof, 1));
            %                 for i = 1:obj.Ndof
            %                     for j = 1:obj.Ndof
            %                         corilisVec3(i) = corilisVec3(i) - dYvar(j)*1/2*gradM(k + (j-1)*obj.Ndof, i)*dYvar(k);
            %                     end
            %                 end
            %                 Ce3{k} = corilisVec3;
            %             end
            %             for j = 1:obj.Ndof
            %                 Ce3sum = Ce3sum + Ce3{k};
            %             end
            %             CeSum = Ce1sum+Ce2sum+Ce3sum;
            %%%%%%%%%%%%%%%%%%% slow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     corilisVec1(i) = corilisVec1(i) + dYvar(j)*1/2* gradM(j + (i-1)*obj.Ndof, :)*dYvar;
            %                 end
            %             end
            %             for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                     corilisVec2(i) = corilisVec2(i) + dYvar(j)*1/2* gradM( (1:obj.Ndof) + (i-1)*obj.Ndof, j)*dYvar;
            %                 end
            %             end
            %              for i = 1:obj.Ndof
            %                 for j = 1:obj.Ndof
            %                      corilisVec3(i) = corilisVec3(i) - dYvar(j)*1/2*gradM((1:obj.Ndof)  + (j-1)*obj.Ndof, i)*dYvar;
            %                 end
            %             end
            %%%%%% sanity check %%%%%%%%%%%%%%%%%
            %   C  = A'*(obj.robot.Mmat*Adot + obj.robot.corilos*A);
            %   J = Jori*A;
            % G = A'*Gori
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            expr = {};
            expr{end+1} = SymFunction('Ce',eval_math_fun('InertiaToCoriolis',{obj.robot.Mmat,obj.q,obj.dq}, [],'DelayedSet',false),{obj.q,obj.dq});
            
            path = [obj.export_path, '/+Normal'];
            expr{end+1} = SymFunction('YnormalofQ', Y, obj.q);
            expr{end+1}  = SymFunction('dYnormalOfX', dY, obj.x);
            
            expr{end+1}  = SymFunction('MnormalofQ', Mnormal, obj.q);
            expr{end+1}  = SymFunction('BnormalofQ', Bnormal, obj.q);
            %   expr{end+1}  = SymFunction('CnormalofQ', CeSum, obj.q, dYvar);
            
            % expr{end+1}  = SymFunction('CmatNormalofQ', corilisMat, obj.q, dYvar);
            expr{end+1}  = SymFunction('Gnormal', Gnormal, Yvar);
            expr{end+1}  = SymFunction('J_leftfoot_normal', J_leftFoot, Yvar);
            expr{end+1}  = SymFunction('J_rightfoot_normal', J_rightFoot, Yvar);
            
            expr{end+1}  = SymFunction('JnormalOfQ', JY, obj.q);
            expr{end+1}  = SymFunction('dJnormalOfX', dJY, obj.x);
            
            expr{end+1}  = SymFunction('AnormalofQ', A, obj.q);
            %  expr{end+1}  = SymFunction('AnormalDotOfX', Adot, obj.x);
            expr{end+1}  = SymFunction('AnormalDotOfX', Adot, obj.x);
            
            obj.processExport(expr, path);
        end
        
        function CeSum = CalculateNormalCoriolis(obj, Mnormal, A)
            Ce1 = cell(obj.Ndof,1);
            Ce2 = cell(obj.Ndof,1);
            Ce3 = cell(obj.Ndof,1);
            Ce1sum = SymExpression(zeros(obj.Ndof, 1));
            Ce2sum = SymExpression(zeros(obj.Ndof, 1));
            Ce3sum = SymExpression(zeros(obj.Ndof, 1));
            
            for k = 1:obj.Ndof
                corilisVec1 = SymExpression(zeros(obj.Ndof, 1));
                for i = 1:obj.Ndof
                    for j = 1:obj.Ndof
                        gradM = jacobian(Mnormal(i,j), obj.q)*A(:,k);
                        corilisVec1(i) = corilisVec1(i) + dYvar(j)*1/2* gradM*dYvar(k);
                    end
                end
                Ce1{k} = corilisVec1;
            end
            
            for k = 1:obj.Ndof
                corilisVec2 = SymExpression(zeros(obj.Ndof, 1));
                for i = 1:obj.Ndof
                    for j = 1:obj.Ndof
                        gradM = jacobian(Mnormal(i,k), obj.q)*A(:,j);
                        corilisVec2(i) = corilisVec2(i) + dYvar(j)*1/2* gradM*dYvar(k);
                    end
                end
                Ce2{k} = corilisVec2;
            end
            
            for k = 1:obj.Ndof
                corilisVec3 = SymExpression(zeros(obj.Ndof, 1));
                for i = 1:obj.Ndof
                    for j = 1:obj.Ndof
                        gradM = jacobian(Mnormal(j,k), obj.q)*A(:,i);
                        corilisVec3(i) = corilisVec3(i) - dYvar(j)*1/2*gradM*dYvar(k);
                    end
                end
                Ce3{k} = corilisVec3;
            end
            
            for j = 1:obj.Ndof
                Ce1sum = Ce1sum + Ce1{k};
                Ce2sum = Ce2sum + Ce2{k};
                Ce3sum = Ce3sum + Ce3{k};
            end
            
            CeSum = Ce1sum + Ce2sum + Ce3sum;
        end
        
        function CeSum = CalculateNormalCoriolisRobotLink(obj, M, A, dY)
            N = obj.Ndof;
            Ce1 = cell( N,1);
            Ce2 = cell( N,1);
            Ce3 = cell( N,1);
            Ce1sum = SymExpression(zeros( N,1));
            Ce2sum = SymExpression(zeros( N,1));
            Ce3sum = SymExpression(zeros( N,1));
            q = obj.q;
            for j = 1:N
                corilisVec1 = SymExpression(zeros(N, 1));
                for i = 1:N
                    for k = 1:N
                        i,j,k
                        corilisVec1(i) = corilisVec1(i) + 1/2*dY(k)*eval_math_fun('Dot', {jacobian(M(i,j), q), A(:,k)});
                    end
                end
                corilisVec1 = corilisVec1*dY(j);
                Ce1{j} = corilisVec1;
            end
            
            for j = 1:N
                corilisVec2 = SymExpression(zeros(N, 1));
                for i = 1:N
                    for k = 1:N
                        i,j,k
                        corilisVec2(i) = corilisVec2(i) + 1/2*dY(k)*eval_math_fun('Dot',{jacobian(M(i,k), q), A(:,j)});
                    end
                end
                corilisVec2 = corilisVec2*dY(j);
                Ce2{j} = corilisVec2;
            end
            
            for j = 1:N
                corilisVec3 = SymExpression(zeros(N, 1));
                for i = 1:N
                    for k = 1:N
                        i,j,k
                        corilisVec3(i) = corilisVec3(i) - 1/2*dY(k)*eval_math_fun('Dot',{jacobian(M(j,k), q), A(:,i)});
                    end
                end
                corilisVec3 = corilisVec3*dY(j);
                Ce3{j} = corilisVec3;
            end
            
            for j = 1:N
                Ce1sum = Ce1sum + Ce1{k};
                Ce2sum = Ce2sum + Ce2{k};
                Ce3sum = Ce3sum + Ce3{k};
            end
            
            CeSum = Ce1sum + Ce2sum + Ce3sum;
        end
        
        function NewFootExpressions(obj)
            RobotModel = obj.robot;
            X = obj.q; dX = obj.dq; Export_path = obj.export_path;
            
            %% foot heel and toes
            Left_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.LeftToe);
            Right_toe = RobotModel.getCartesianPosition(RobotModel.ContactPoints.RightToe);
            
            J_Left_toe = jacobian(Left_toe, obj.q);
            J_Right_toe = jacobian(Right_toe, obj.q);
            
            dJ_Left_toe = SymExpression(zeros(3, obj.Ndof));
            dJ_Right_toe = SymExpression(zeros(3, obj.Ndof));
            
            for i = 1:3
                for j = 1: obj.Ndof
                    dJ_Left_toe(i,j) = jacobian(J_Left_toe(i,j),obj.q)*obj.dq;
                    dJ_Right_toe(i,j) = jacobian(J_Right_toe(i,j),obj.q)*obj.dq;
                end
            end
            
            expr = {};
            expr{end+1}  = SymFunction('dJ_leftToeNew',dJ_Left_toe,obj.x);
            expr{end+1} = SymFunction('dJ_rightToeNew',dJ_Right_toe,obj.x);
            
            path = [obj.export_path, '/+Normal'];
            obj.processExport(expr, path);
        end
        
        
        
    end
    
    methods
        function processExport(obj, expr, path)
            % Export the files
            if ~exist(path,'dir')
                mkdir(char(path));
            end
            
            for i = 1:length(expr)
                expr{i}.export(path);
            end
        end
    end
    
end

