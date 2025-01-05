classdef EOM_walker < handle
    %EOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Mass
        Bu
        Gvec
        h % h = Coriolis + Gravity 
        CGvec % sum of Coriolis and Gvec
        B_hat  % 
        fx_hat_bottom % bottom part from Equation of Motion Dddq + h = JT F + Bu 
        gx
        fx
        gx_hat
        gx_hat_bottom
        
        J_hol_leftContact % include rods and left foot contact
        dJ_hol_leftContact
        J_hol_rightContact
        dJ_hol_rightContact
        
        J_leftStance
        dJ_leftStance
        
        J_rightStance
        dJ_rightStance
        
        pToeHeel % was used for spring-damper model of ground
        Jtoeheel  % may be useful for simulation of continuous controller
        vToeHeel
        dJtoeheel
        aToeHeel
        %
        J_hol_domain
        dJ_hol_domain
        J_swing 
        dJ_swing
        q
        dq

        QPvariableType % 'onlyU' 'U-Fhol', 'U-Fhol-ddq'
        
        A_Fhol2U  = []
        b_Fhol2U  = []
        
        Contact
        ndof = 7
        
        invM
        xinv
        partial_fx_hat_bottom_q
        partial_gx_hat_bottom_q;
        partial_fx_hat_bottom_dq
        partial_gx_hat_bottom_dq;
            
        SVA
    end
    
    properties
        motorIdx = 4:7;
        leftHip = 6; 
        leftKnee = 7; 
        rightHip = 4; 
        rightKnee = 5;
    end
    
    methods
        
        function obj = EOM_walker(q_size, eomtype)
            % varargin for spring stiffness and dampings
            obj.Mass = zeros(q_size);
            obj.invM = zeros(q_size); 
            obj.Bu = u_map_five_link_walker(zeros(q_size,1));
            obj.Gvec = zeros(q_size,1);
            obj.CGvec = zeros(q_size,1);
            obj.B_hat = obj.Bu;
            obj.QPvariableType = eomtype.QPvariableType;  
            
            temp =load('fiveLinkWalkerOriSVA.mat');
            obj.SVA = temp.fiveLinkWalkerSVAOri;
            
%              temp =load('fiveLinkWalkerSVA.mat');
%             obj.SVA = temp.fiveLinkWalkerSVA;
        end

        function obj = updateWalking(obj, q, dq, Contact)
            % defalut using midfoot as contact 
            obj.q = q;
            obj.dq = dq;
            obj.updateEOM(q,dq); 
            obj.updateHolJacobian(q,dq);
            
            x = [q;dq];
            % augmented dynamics
            obj.h = obj.CGvec; 
            obj.Contact = Contact; 
            if strcmp(Contact,      'LeftFootContact')
                obj.B_hat = [ obj.Bu, obj.J_leftStance', zeros(size(obj.J_rightStance'))];
                obj.J_hol_domain = [obj.J_leftStance];
                obj.dJ_hol_domain = [ obj.dJ_leftStance];    
                obj.J_swing = obj.J_rightStance; 
                obj.dJ_swing = obj.dJ_rightStance;
            elseif strcmp(Contact,  'RightFootContact')
                obj.B_hat = [ obj.Bu, zeros(size(obj.J_leftStance')),obj.J_rightStance'];
                obj.J_hol_domain =  [ obj.J_rightStance];
                obj.dJ_hol_domain = [ obj.dJ_rightStance];    
                obj.J_swing = obj.J_leftStance; 
                obj.dJ_swing = obj.dJ_leftStance;
            end
            obj.invM = inv(obj.Mass); 
            obj.fx_hat_bottom = obj.invM*(-obj.CGvec);
            obj.gx_hat_bottom = obj.invM*obj.B_hat;
            obj.fx = [dq; obj.fx_hat_bottom];
            
            if strcmp(obj.QPvariableType, 'onlyU') || strcmp(obj.QPvariableType, 'OSC')
                %%%%%%%%%%  QP of using only torque as optimization variable
                obj.xinv = inv(obj.J_hol_domain*(obj.invM*obj.J_hol_domain'));
                 %%%%%%%% Fhol = A_Fhol2U*u + b_Fhol2U
                obj.A_Fhol2U = -obj.xinv*obj.J_hol_domain*(obj.invM*obj.Bu);
                obj.b_Fhol2U = obj.xinv*(obj.J_hol_domain*(obj.invM*obj.h) - obj.dJ_hol_domain*obj.dq);
                obj.B_hat = obj.Bu + obj.J_hol_domain'*obj.A_Fhol2U; 
                obj.fx_hat_bottom = obj.fx_hat_bottom + obj.invM*(obj.J_hol_domain'*obj.b_Fhol2U); 
                obj.gx_hat_bottom = obj.invM*obj.B_hat;  
            end
        end
    end
    
    
    methods %%%%% gradients
        %%%%%% note that the Jacobain is different for calculating the
        %         %%%%%% gradients of the vector fields vs. that of the reset map
        
        %%%%% first calculate the gradient of the vector field fx and gx.
        
        function cal_partial_dynamics(obj, q, dq, stanceLeg)
            
            %%%% require
            gradM = obj.evaluateGradM(q, dq);
            
            if strcmp(stanceLeg, 'RightFootContact')
                [grad_J_q, grad_J_dq] = obj.evaluateGrad_J_right(q,dq); %%% left swing
                [grad_dJ_q, grad_dJ_dq] = obj.evaluateGrad_dJ_right(q,dq); %%% left swing
            else 
                  [grad_J_q, grad_J_dq] = obj.evaluateGrad_J_left(q,dq); %%% left swing
                [grad_dJ_q, grad_dJ_dq] = obj.evaluateGrad_dJ_left(q,dq); %%% left swing
            end
            
            gradInvM  = -obj.MmT(obj.invM, obj.TmM(gradM, obj.invM));
            J = obj.J_hol_domain;
            gradxinv = obj.gradxinv_q(J, gradInvM,  grad_J_q); %2.2.n
            
            grad_h_q = grad_CG_q([q;dq]);
            grad_h_dq = grad_CG_dq([q;dq]);
            
            grad_A_Fhol2U_q = obj.grad_A_Fhol2u( gradxinv, grad_J_q,gradInvM); %2.4.n
            [grad_b_Fhol2U_q , grad_b_Fhol2U_dq] = obj.grad_b_Fhol2U(dq, gradxinv, gradInvM, grad_h_q,grad_h_dq, grad_dJ_q, grad_dJ_dq);
            
            obj.partial_fx_hat_bottom_q = obj.gradfx_q( gradM, grad_J_q, grad_b_Fhol2U_q, grad_h_q);
            obj.partial_fx_hat_bottom_dq = obj.gradfx_dq(grad_b_Fhol2U_dq, grad_h_dq );
            [obj.partial_gx_hat_bottom_q, obj.partial_gx_hat_bottom_dq] = ...
                obj.gradgx_x( gradInvM, grad_J_q, grad_A_Fhol2U_q);
        end
        
        function out = gradfx_q(obj, gradM, gradJ, grad_b_Fhol2U_q, grad_h_q) 
            %%% gradM n.n.n 
            %%% gradJ 2.n.n 
            %%% grad_b_Fhol2U_q 2.n
            %%% out = n.n
            %%%% J is the stance Jacobain %%%%%%%%%%%%%
            out = squeeze(obj.TmM(gradM, (-obj.CGvec))) + obj.invM*(-grad_h_q); 
            if strcmp(obj.QPvariableType, 'onlyU') || strcmp(obj.QPvariableType, 'OSC')
                out = out + squeeze(obj.TmM(gradM, (obj.J_hol_domain'*obj.b_Fhol2U))) + ...
                        obj.invM*( squeeze(obj.TmM(permute(gradJ, [2,1, 3]), obj.b_Fhol2U)) + ... %%% check permute
                        obj.J_hol_domain'*grad_b_Fhol2U_q);
            end
        end 
        
        function out = gradfx_dq(obj, grad_b_Fhol2U_dq, grad_h_dq )
            %%% grad_b_Fhol2U_dq : 2.n
            %%%% out = n.n
            %%%
            out = obj.invM*(-grad_h_dq); 
            if strcmp(obj.QPvariableType, 'onlyU') || strcmp(obj.QPvariableType, 'OSC')
                out = out +  obj.invM*(obj.J_hol_domain'*grad_b_Fhol2U_dq);
            end
        end
        
        function [out_q, out_dq] = gradgx_x(obj, gradInvM, grad_J_q, grad_A_Fhol2U_q)
            %%% gradInvM : n.n.n
            %%% grad_J_q 2.n.n
            %%% obj.A_Fhol2U: 2*4
            %%% obj.Bhat: n.4 (4 inputs)
            %%% grad_A_Fhol2U_q: 2*4*n
            %%% out_q: n.4.n
            out_q =obj.TmM(gradInvM, obj.B_hat) + ...
                obj.MmT(obj.invM, obj.TmM(permute(grad_J_q, [2,1, 3]), obj.A_Fhol2U) ) + ...
                + obj.MmT(obj.invM*obj.J_hol_domain', grad_A_Fhol2U_q);
            out_dq = zeros(size(out_q)); %obj.invM*obj.J_hol_domain'*grad_A_Fhol2U_dq;
        end
         
        function [out_q] = grad_A_Fhol2u(obj, gradxinv, grad_J_q,GradInvM) 
            %%% xinv. 2.2 
            %%%gradxinv 2.2.n
            %%% obj.Bu n.4
            %%%grad_J_q: 2.n.n
            %%% out_q 2.4.n
            out_q = obj.TmM(gradxinv, obj.J_hol_domain*(obj.invM*obj.Bu)) + ...
                      obj.MmT(obj.xinv, obj.TmM(grad_J_q, obj.invM*obj.Bu)) + ...
                      obj.MmT((obj.xinv*obj.J_hol_domain),obj.TmM(GradInvM, obj.Bu));
           % grad_dq = zeros(size(grad_q)); 
        end 
                 
        function [out_q, out_dq] = grad_b_Fhol2U(obj, dq, gradxinv, gradInvM, grad_h_q,grad_h_dq, grad_dJ_q, grad_dJ_dq)
            %%% grad_dJ_q, grad_dJ_dq: 2.n.n
            %%% grad_h_q, grad_h_dq: n.n
            %%% gradxinv %2.2.n
            %%% out: 2.n
           out_q = squeeze(obj.TmM(gradxinv, (obj.J_hol_domain*(obj.invM*obj.h)) - obj.dJ_hol_domain*obj.dq))+ ...
               obj.xinv*squeeze(obj.TmM(grad_dJ_q, obj.invM*obj.h)) + ...
               obj.xinv*obj.J_hol_domain*squeeze(obj.TmM(gradInvM, obj.h)) + ...
               obj.xinv*obj.J_hol_domain*obj.invM*grad_h_q -  obj.xinv*squeeze(obj.TmM(grad_dJ_q,dq));
            
           out_dq =  obj.xinv*obj.J_hol_domain*(obj.invM*grad_h_dq) - ...
                        obj.xinv*obj.dJ_hol_domain - obj.xinv*squeeze(obj.TmM(grad_dJ_dq, dq));
        end
         
        function  out = gradxinv_q(obj, J, gradInvM,  grad_J_q)
            %%%%% calculate grad(invM, q), grad( inv(J*inv(M)*J'), q)
            invM = obj.invM;
            A2 = inv(J*invM*J');
            %%% gradJMJ: zeros(2, 2, obj.ndof);
            %%% gradxinv = zeros(2, 2, obj.ndof);
            gradJMJ = obj.TmM(grad_J_q, invM*J') + obj.MmT(J, obj.TmM(gradInvM,J'))...
                + obj.MmT(J*invM, permute(grad_J_q, [2,1, 3]));
            out = -obj.MmT(A2, obj.TmM(gradJMJ, A2));
        end
%         
        function [gradM] = evaluateGradM(obj, q,dq)
            %gradM = reshape(grad_Mq(q), );
            %%% to avoid erros, manully assign
            gradM = zeros(obj.ndof, obj.ndof, obj.ndof);
            stackedGradM = grad_Mq([q;dq]);
            for i = 1:obj.ndof
                for j = 1:obj.ndof
                    gradM(i, j, :) = stackedGradM(j+ (i-1)*obj.ndof,:);
                end
            end
        end
%         
        function [ grad_J_q, grad_J_dq] = evaluateGrad_J_left(obj, q,dq)
            %%%% which Jacobian is this.
            %gradJ = reshape(grad_Jq(q), );
            %%% to avoid erros, manully assign
            grad_J_q = zeros(2, obj.ndof, obj.ndof);
            grad_J_dq = zeros(2, obj.ndof, obj.ndof);
            
            stackedGradq = grad_J_left_q([q;dq]);
            
            for i = 1:2
                for j = 1:obj.ndof
                    grad_J_q(i, j, :) = stackedGradq(j+ (i-1)*obj.ndof,:);
                end
            end
        end
        
        function [grad_dJ_q, grad_dJ_dq] = evaluateGrad_dJ_left(obj, q,dq)
            %%%% which Jacobian is this.
            %gradJ = reshape(grad_Jq(q), );
            %%% to avoid erros, manully assign
            grad_dJ_q = zeros(2, obj.ndof, obj.ndof);
            grad_dJ_dq = zeros(2, obj.ndof, obj.ndof);
            
            stackedGradq = grad_dJ_left_q([q;dq]);
            stackedGraddq = grad_dJ_left_dq([q;dq]);
            
            for i = 1:2
                for j = 1:obj.ndof
                    grad_dJ_q(i, j, :) = stackedGradq(j+ (i-1)*obj.ndof,:);
                end
            end
            for i = 1:2
                for j = 1:obj.ndof
                    grad_dJ_dq(i, j, :) = stackedGraddq(j+ (i-1)*obj.ndof,:);
                end
            end
        end
        
        function [ grad_J_q, grad_J_dq] = evaluateGrad_J_right(obj, q,dq)
            %%%% which Jacobian is this.
            %gradJ = reshape(grad_Jq(q), );
            %%% to avoid erros, manully assign
            grad_J_q = zeros(2, obj.ndof, obj.ndof);
            grad_J_dq = zeros(2, obj.ndof, obj.ndof);
            
            stackedGradq = grad_J_right_q([q;dq]);
            %stackedGraddq = grad_J_right_dq([q;dq]);
            
            for i = 1:2
                for j = 1:obj.ndof
                    grad_J_q(i, j, :) = stackedGradq(j + (i-1)*obj.ndof,:);
                end
            end
        end
        
        function [ grad_dJ_q, grad_dJ_dq] = evaluateGrad_dJ_right(obj, q,dq)
            %%%% which Jacobian is this.
            %gradJ = reshape(grad_Jq(q), );
            %%% to avoid erros, manully assign
            grad_dJ_q = zeros(2, obj.ndof, obj.ndof);
            grad_dJ_dq = zeros(2, obj.ndof, obj.ndof);
            
            stackedGradq = grad_dJ_right_q([q;dq]);
            stackedGraddq = grad_dJ_right_dq([q;dq]);
            
            for i = 1:2
                for j = 1:obj.ndof
                    grad_dJ_q(i, j, :) = stackedGradq(j+ (i-1)*obj.ndof,:);
                end
            end
            for i = 1:2
                for j = 1:obj.ndof
                    grad_dJ_dq(i, j, :) = stackedGraddq(j+ (i-1)*obj.ndof,:);
                end
            end
        end
          
        function out = gx_hat_bottom_func(obj, x) 
            %%% for some fake state x
            q = x(1:obj.ndof); 
            dq = x(obj.ndof+1:end);
            Mass = Mmat_five_link_walker(q);
            %%% right foot stance
            J_hol = J_rightToe(q);
            J_hol = [J_hol(1,:); J_hol(2,:)]; 
            xinv = inv(J_hol*(Mass\J_hol'));
            %%%%%%%% Fhol = A_Fhol2U*u + b_Fhol2U
            A_Fhol2U = -xinv*J_hol*(Mass\obj.Bu);
            B_hat = obj.Bu + J_hol'*A_Fhol2U;
            out = Mass\B_hat;
        end 
    end
    
    methods %(Access = private)
        function Coriolis = coriolis(obj, q, dq) 
            q_size = 7;
            C = zeros(q_size,1);
            %%% any speed-up by commenting this ?
            for i = 1:3
                for j = 1:q_size
                    C = C + feval(strcat('Ce',int2str(i),'_vec',int2str(j),'_five_link_walker'),q,dq);  %~9ms;
                end
            end
            Coriolis = -C;
        end 
        
        function updateEOM(obj,q,dq)
            obj.Mass = Mmat1_Robot_Assembly_v3_straight_leg(q)+Mmat2_Robot_Assembly_v3_straight_leg(q)+Mmat3_Robot_Assembly_v3_straight_leg(q)...
                    + Mmat4_Robot_Assembly_v3_straight_leg(q)+ Mmat5_Robot_Assembly_v3_straight_leg(q) +Mmat6_Robot_Assembly_v3_straight_leg(q);
            obj.CGvec = -CGvec_five_link_walker(q, dq);
        end
        
        function updateHolJacobian(obj,q,dq)
            %% update Jacobians
            % foot jacobians
            obj.J_rightStance = J_rightToe(q);
            obj.J_leftStance =  J_leftToe(q);
            obj.dJ_rightStance = dJ_rightToe([q;dq]);
            obj.dJ_leftStance =  dJ_leftToe([q;dq]);
        end
        
    end
 
    methods %%% tensor matrix multiplcation
        
        function out = TmM(obj, T, M)
            %%% tensor x matrix 
            [t1, t2, t3] = size(T); 
            [m1, m2] = size(M); 
            %%% t2 should be the same as m1
            out = zeros(t1, m2, t3); 
            if t2 == m1
                for i = 1:t3
                    out(:, :, i) = T(:,:, i)*M;
                end
            else 
                disp('wrong dimension');
            end
        end
        
        function out = MmT(obj, M, T)
            %%% matrix x tensor
            [m1, m2] = size(M);
            [t1, t2, t3] = size(T); 
            %%% m2 should be the same as t1
            out = zeros(m1, t2, t3); 
            if m2 == t1
                for i = 1:t3
                    out(:, :, i) = M*T(:,:, i);
                end
            else 
                disp('wrong dimension');
            end
        end 
        
    end 
    
end