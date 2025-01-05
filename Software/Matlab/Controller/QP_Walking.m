classdef QP_Walking < handle

    properties
        %% QP parameters
        n_u
        n_GRF
        n_hol
        n_rod 
        n_delta
        n_u_hat
        n_total   % total variables in qp
        %%
        % setting for ZMP stability
        AneqGRF
        BneqGRF
        AneqGRFleft
        BneqGRFleft
        AneqGRFright
        BneqGRFright
       
        % this should be only used for heelToe two point contact
        AeqGRFx
        BeqGRFx  % toe heel x force should be identical as foot is non-compressible
        
        AeqGRF %% setting for force control
        BeqGRF
        AeqGRFleft
        BeqGRFleft
        AeqGRFright
        BeqGRFright
        % properities for formulating QP
        Aeq
        Beq
        Aneq
        Bneq
        H
        g
        lb
        ub
        
        Ahol 
        Bhol
        % different QP formulation
        Beq_eom
        Aeq_eom
        
        % holonomic rod
        Beq_rod
        Aeq_rod
        
        Beq_2D
        Aeq_2D
        % left contact hol
        Beq_left
        Aeq_left
        % right contact hol
        Beq_right
        Aeq_right
        Aclf
        Bclf
            
        %% external parameters
        motorTorqueLimits = 2000*ones(1,4)  
        miu
         
        ddq_size = 7
        contact 
        initialGuess 
        ContactMatrix
        QPvariableType %%% 'onlyU' 'U-Fhol', 'U-Fhol-ddq'
    end
    
    methods
        function obj = QP_Walking(miu, QPvariableType)
            obj.miu = miu;
            % actuator limites - u_max - u <=0;  - u_max + u <=0;
            obj.lb = -obj.motorTorqueLimits';
            obj.ub = -obj.lb;
            obj.initialGuess = [];
            obj.QPvariableType = QPvariableType;
            obj.ConfigMidFootQP()
        end
        
        function exitflag = convertMatlabQuadprog2QPOASES(obj,exitflag)
            % just to make the exitflag consistent change from Quadprog to
            % QPOASES
            switch exitflag
                case 1
                    exitflag = 0; % solved
                case 0
                    exitflag = 1; % num of iteration exceeded
                case -2
                    exitflag = -2; % infeasible
                case -3
                    exitflag = -3; % unbounded
                otherwise
                    exitflag = exitflag;
            end
        end
        
        function [A,lbA, ubA] = convertConstraints(obj)
            A = [obj.Aneq; obj.Aeq];
            ubA = [obj.Bneq; obj.Beq];
            lbA = [-inf*ones(size(obj.Bneq));obj.Beq];
        end
        
        function ConfigMidFootQP(obj)
            switch obj.QPvariableType %'onlyU' 'U-Fhol', 'U-Fhol-ddq'
                case 'onlyU'
                    obj.ConfigMidFootQP_onlyU();
                case 'U-Fhol'
                    obj.ConfigMidFootQP_withFhol();
                case 'U-Fhol-ddq'
                    obj.ConfigMidFootQP_withddq();
                case 'OSC'
                    obj.ConfigMidFootQP_OSC();
                otherwise
                    disp('wrong QP type')
            end
        end
        
        function [u,exitflag,F_GRF, delta, fval, numiter, QPfail] = constructAndSolve(obj, ...
                   contact ,uf0, eom, outputs)
            % put the safesolve here
            switch obj.QPvariableType %'onlyU' 'U-Fhol', 'U-Fhol-ddq'
                case 'onlyU'
                    obj.configureCost_onlyU( A_hat, QP_Lf, smoothingQP,contact, eom);
                    obj.ConfigureConstraints_onlyU( A_hat,QP_Lf,psi0,psi1,eom, contact);
                    [u,exitflag,F_GRF,delta, fval, numiter] = obj.solve_onlyU(eom);
                case 'U-Fhol'
                    obj.configureCost_withFhol( A_hat, QP_Lf, smoothingQP,contact);
                    obj.ConfigureConstraints_withFhol( A_hat, QP_Lf, psi0, psi1, eom, contact);
                    [u,exitflag,F_GRF,delta, fval, numiter] = obj.solve_withFhol();
                case 'U-Fhol-ddq'
                    obj.configureCost_withddq( A_hat, QP_Lf, smoothingQP,contact);
                    obj.ConfigureConstraints_withddq( A_hat, QP_Lf, psi0, psi1, eom, contact);
                    [u, exitflag, F_GRF, ~, delta, fval, numiter] = obj.solve_withddq();
                case 'OSC'
                    obj.configureCost_OSC(outputs, eom);
                    obj.ConfigureConstraints_OSC(eom, contact);
                    [u,exitflag,F_GRF, fval, numiter] = obj.solve_OSC(eom);
                    delta = 0; 
                otherwise
                    disp('wrong QP type')
            end
            
            if exitflag ~= 0 %&& d == 1
                switch obj.QPvariableType %'onlyU' 'U-Fhol', 'U-Fhol-ddq'
                    case 'onlyU'
                        [u,exitflag, F_GRF, delta, fval, numiter] = obj.safesolve_onlyU(uf0, eom);
                    case 'U-Fhol'
                        [u,exitflag, F_GRF, delta, fval, numiter] = obj.safesolve_withFhol(uf0);
                    case 'U-Fhol-ddq'
                        [u, exitflag, F_GRF, ~, delta, fval, numiter] = obj.safesolve_withddq(uf0);
                    case 'OSC'
                        [u, exitflag, F_GRF, delta, fval,numiter] =obj.safesolve_OSC(uf0, eom);
                    otherwise
                        disp('wrong QP type')
                end
                exitflag = exitflag + 10;
                if exitflag~= 10
                    u = zeros(4,1);
                end
                QPfail = 1;
            else
                QPfail = 0;
            end
        end
        
        function [u,exitflag,Frod,F_GRF, delta, fval, numiter, QPfail] = feedbackLinearization(obj, outputs, A_hat, Lf, eom, contact) 
            % should be quick and easy as a continuous closed-form feedback
            % law % mainly for debugging purposes
            RD1 = outputs.RD1;
            RD2 = outputs.RD2;
            if strcmp(obj.QPvariableType, 'onlyU') %  only for this formulation 
                % 1Dhopping
                eplison = 50; 
                switch contact
                    case 'Prejump' % 1 RD1 output, 5 RD2 outputs
                        v = [-eplison*RD1.y; 
                             -2*eplison*RD2.dy - eplison^2*RD2.y];
                    case 'Flight' % 10 RD2 outputs
                          v = -2*eplison*RD2.dy - eplison^2*RD2.y;
                    case 'Landing' % 6 RD2 outputs]
                          v = -2*eplison*RD2.dy - eplison^2*RD2.y;
                end
                u = A_hat\(-Lf + v);

                Fhol = eom.A_Fhol2U*u + eom.b_Fhol2U;
                Frod = Fhol(1:2); 
                if length(Fhol) >2 % on the ground
                    F_GRF = Fhol(3:end); 
                else
                    F_GRF = zeros(10,1); % off the ground
                end
                % just for consistency of the code
                exitflag = 0; delta = 0; fval= 0;  numiter = 0;  QPfail = 0;
            end
        end
        
    end
    
    methods % conventional CLF-QP formulation for robots
        
        function ConfigMidFootQP_withFhol(obj)
            obj.n_u = 10; obj.n_rod = 2; 
            obj.n_GRF = 10;
            obj.n_delta = 1; 
            if strcmp(obj.model2D3D, '2D')
                obj.hol_size_2D = 3;
            else
                obj.hol_size_2D = 0;
            end
 
            n_arguments = [obj.n_u, obj.n_rod, obj.n_GRF, obj.hol_size_2D, obj.n_delta];
            obj.n_total = sum(n_arguments); 
            obj.n_u_hat = obj.n_total - 1; 
            
            obj.lb = [-obj.motorTorqueLimits;-inf*ones(obj.n_total-obj.n_u,1)];
            obj.ub = -obj.lb;
            obj.Aeq = [];
            obj.Beq = [];
            obj.Aneq = [];
            obj.Bneq = [];
            % Equality constraints on GRF: identity matrix for doing force control (desired Force)
            obj.AeqGRFleft = [zeros(obj.n_GRF/2,obj.n_u + obj.n_rod),eye(obj.n_GRF/2),...
                zeros(obj.n_GRF/2,obj.n_GRF/2+ obj.hol_size_2D + obj.n_delta)];
            obj.BeqGRFleft = zeros(obj.n_GRF/2,1);
            obj.AeqGRFright = [zeros(obj.n_GRF/2,obj.n_u + obj.n_rod + obj.n_GRF/2),eye(obj.n_GRF/2),...
                zeros(obj.n_GRF/2,obj.hol_size_2D + obj.n_delta)];
            obj.BeqGRFright = zeros(obj.n_GRF/2,1);
            
            % friction cone u_friction_cone_RightMidFoot(RightForceXYZ, miu)>=0;
            % moment constraint; not consider yaw moment for now.
            %appx_miu = sqrt(2)/2*obj.miu; % polyhyd appx
            appx_miu = obj.miu; % 
            
            if strcmp(obj.model2D3D ,'3D')
                obj.ContactMatrix = [
                    1 0 -appx_miu 0 0;
                    -1 0 -appx_miu 0 0;   % forward friction
                    0 1 -appx_miu 0 0;
                    0 -1 -appx_miu 0 0;   % lateral friction
                    0 0 -1 0 0;           % normal force
                    0 0 -obj.half_foot_length 1 0;
                    0 0 -obj.half_foot_length -1 0;  %forward zmp
                    0 0 -obj.half_foot_length*obj.miu 0 1;
                    0 0 -obj.half_foot_length*obj.miu 0 -1;
                    ];
            elseif strcmp(obj.model2D3D ,'2D')
                obj.ContactMatrix = [
                 %   1 0 -appx_miu 0 0; % due to the holonomic constraint   on the pelvis, this friction has to be commented out
                 %  -1 0 -appx_miu 0 0;   % forward friction % this needs to be removed due to  
                    0 0 -1 0 0;           % normal force
                  0 0 -obj.half_foot_length 1 0;
                   0 0 -obj.half_foot_length -1 0;  %forward zmp
                    ];
            elseif strcmp(obj.model2D3D ,'2Dsoft')
                obj.ContactMatrix = [
             %       1 0 -appx_miu 0 0;
               %    -1 0 -appx_miu 0 0;   % forward friction
                    0 0 -1 0 0;           % normal force
                    0 0 -obj.half_foot_length 1 0;
                    0 0 -obj.half_foot_length -1 0;  %forward zmp
                    ];
            end
            rowsGRF = size(obj.ContactMatrix,1);
            
            % GRF force Fznormal>0 frictionxy<uFnorm, My <lt lh*Fz Mz <lt lh * Fy
            obj.AneqGRF = [zeros(rowsGRF,obj.n_u + obj.n_rod), obj.ContactMatrix,zeros(rowsGRF,obj.n_GRF/2 + obj.hol_size_2D + obj.n_delta);
                           zeros(rowsGRF,obj.n_u + obj.n_rod + obj.n_GRF/2), obj.ContactMatrix,zeros(rowsGRF, obj.hol_size_2D + obj.n_delta)];
            obj.BneqGRF = zeros(2*rowsGRF,1);
            
            obj.AneqGRFleft = [zeros(rowsGRF,obj.n_u+obj.n_rod),obj.ContactMatrix,zeros(rowsGRF,obj.n_GRF/2 + obj.hol_size_2D + obj.n_delta)];
            obj.BneqGRFleft = zeros(rowsGRF,1);
            obj.AneqGRFright = [zeros(rowsGRF,obj.n_u +obj.n_rod + obj.n_GRF/2),obj.ContactMatrix,zeros(rowsGRF,obj.hol_size_2D + obj.n_delta)];
            obj.BneqGRFright = zeros(rowsGRF,1);
        end
        
        function configureCost_withFhol(obj, A_hat, QP_Lf, smoothingQP,contact,d)
            ATA = A_hat'*A_hat;
            obj.H = blkdiag(ATA, obj.p);
            obj.g = [(QP_Lf'*A_hat)';zeros(obj.n_delta, 1)];
            
            grfMin = blkdiag(1,1,0,1,1);
            grfMIN = blkdiag(grfMin, grfMin);
            obj.H =  obj.H + smoothingQP*blkdiag( eye(obj.n_u), zeros(obj.n_rod), 0*grfMIN, zeros(obj.hol_size_2D), 0); %smoother and faster
            %obj.H =  blkdiag(eye(obj.n_u),zeros(obj.n_hol),eye(obj.n_GRF),obj.p);
            % if sum(desiredGRF)~= 0
            %obj.H = obj.H + blkdiag(zeros(obj.n_u),obj.P_GRF,zeros(obj.n_hol+1));
            % obj.g = obj.g - [zeros(obj.n_u,1);obj.P_GRF*desiredGRF;zeros(obj.n_hol+1,1)];
            % end
%             obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_hol),0.*eye(obj.n_GRF),obj.p);
%             obj.g = 0.*obj.g + zeros(obj.n_u_hat+1,1);
            % if minimize u; this could work
            if ( strcmp(contact,'DoubleSupportFromLeft') || strcmp(contact,'DoubleSupportFromRight'))  %&& d> 1
                if d == 1
                    obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_rod),0.*eye(obj.n_GRF),0.*eye(obj.hol_size_2D),obj.p);
                    obj.g = 0.*obj.g + zeros(obj.n_u_hat+1,1);
                else
                    obj.H = .1*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_rod),0.*eye(obj.n_GRF),0.*eye(obj.hol_size_2D),obj.p);
                    obj.g = .1*obj.g + zeros(obj.n_u_hat+1,1);
                end
            end
%             obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_rod),0.*eye(obj.n_GRF),0.*eye(obj.hol_size_2D),obj.p);
%             obj.g = 0.*obj.g + zeros(obj.n_u_hat+1,1);
        end
        
        function ConfigureConstraints_withFhol(obj, A_hat, QP_Lf, psi0, psi1, eom, contact)
            obj.contact  = contact; 
            %% equality constraints
            % dJ_hol*dq + J_hol*eom.Mass\(eom.B_hat*u_hat - eom.Coriolis - eom.Gvec + eom.springVec) = 0;
            % this implictly enforce same x force for toe and heel
            Ahol = eom.J_hol_domain*(eom.Mass\eom.B_hat);
            
            obj.Bhol = eom.J_hol_domain*(eom.Mass\( eom.CGvec - eom.springVec)) - eom.dJ_hol_domain*eom.dq;
            if strcmp(obj.model2D3D, '2D') % convert F_additional from u_hat using the additional holonomic constrait relation
                xinv2D = inv(eom.J_2Dpelvis*(eom.Mass\(eom.J_2Dpelvis')));
                Ahol = eom.J_hol_domain*(eom.Mass\eom.B_hat + eom.Mass\(eom.J_2Dpelvis'*xinv2D*eom.J_2Dpelvis*(eom.Mass\eom.B_hat)));
                eomH = eom.CGvec - eom.springVec;
                obj.Bhol = eom.J_hol_domain*(eom.Mass\(eomH - eom.J_2Dpelvis'*xinv2D*eom.J_2Dpelvis*(eom.Mass\eomH))) - eom.dJ_hol_domain*eom.dq;
            end
            obj.Ahol = [Ahol, zeros(length(obj.Bhol),obj.n_delta)];
            
            %% Inequality
            % CLF: psi0 + psi1'*(A_hat*u_hat + Lf) <= delta
            obj.Aclf = [psi1'*A_hat, -1];
            obj.Bclf = -psi1'*QP_Lf - psi0;
            
            switch contact
                case AmberConstants.RightFootContact
                    obj.Aneq = [obj.Aclf;obj.AneqGRFright];
                    obj.Bneq = [obj.Bclf;obj.BneqGRFright];
                    obj.Aeq = [obj.Ahol;  obj.AeqGRFleft];
                    obj.Beq = [obj.Bhol;  obj.BeqGRFleft];
                case AmberConstants.LeftFootContact
                    obj.Aneq = [obj.Aclf;obj.AneqGRFleft];
                    obj.Bneq = [obj.Bclf;obj.BneqGRFleft];
                    obj.Aeq = [obj.Ahol;  obj.AeqGRFright];
                    obj.Beq = [obj.Bhol;  obj.BeqGRFright];
                case 'Flight'
                    obj.Aneq = obj.Aclf ;
                    obj.Bneq = obj.Bclf ;
                    obj.Aeq = [obj.Ahol; obj.AeqGRFleft; obj.AeqGRFright];
                    obj.Beq = [obj.Bhol; obj.BeqGRFleft; obj.BeqGRFright];
                otherwise % 'DoubleSupportFromRight' 'DoubleSupportFromLeft' 'Prejump' 'PeriodicStance' % DSP is the same
                    obj.Aneq = [obj.Aclf;obj.AneqGRF];
                    obj.Bneq = [obj.Bclf;obj.BneqGRF];
                    obj.Aeq = obj.Ahol;
                    obj.Beq = obj.Bhol;
            end
        end
 
        function RemoveConstraints_withFhol(obj)
            switch obj.contact
                case 'DoubleSupportFromRight' % next domain is left foot contact
                    obj.Aneq = [obj.Aclf;obj.AneqGRFleft];
                    obj.Bneq = [obj.Bclf;obj.BneqGRFleft];
                    obj.Aeq = [obj.Ahol;  obj.AeqGRFright];
                    obj.Beq = [obj.Bhol;  obj.BeqGRFright];
                case 'DoubleSupportFromLeft'
                    obj.Aneq = [obj.Aclf;obj.AneqGRFright];
                    obj.Bneq = [obj.Bclf;obj.BneqGRFright];
                    obj.Aeq = [obj.Ahol;  obj.AeqGRFleft];
                    obj.Beq = [obj.Bhol;  obj.BeqGRFleft];
                case 'Prejump'
                    obj.Aneq = obj.Aclf ;
                    obj.Bneq = obj.Bclf ;
                    obj.Aeq = [obj.Ahol; obj.AeqGRFleft; obj.AeqGRFright];
                    obj.Beq = [obj.Bhol; obj.BeqGRFleft; obj.BeqGRFright];
                case 'PeriodicStance'
                    obj.Aneq = obj.Aclf ;
                    obj.Bneq = obj.Bclf ;
                    obj.Aeq = [obj.Ahol; obj.AeqGRFleft; obj.AeqGRFright];
                    obj.Beq = [obj.Bhol; obj.BeqGRFleft; obj.BeqGRFright];
                otherwise 
                    disp('QP fails in SSP')
            end
        end
        
        function [u,exitflag,Frod,F_GRF,delta, fval, numiter] = solve_withFhol(obj)
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 10000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g,A,obj.lb,obj.ub,lbA,ubA,options);
%             options = optimoptions('quadprog','Display','off');
%             [xsol,fval,exitflag, numiter] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,[],options);
            u = xsol(1:obj.n_u);
            Frod = xsol(obj.n_u+1:obj.n_u+obj.n_rod);
            F_GRF = xsol(obj.n_u+obj.n_rod+1:obj.n_u+obj.n_rod+obj.n_GRF);
            delta = xsol(end);
        end
        
        function [u,exitflag,Frod,F_GRF,delta, fval, numiter] = safesolveOld(obj,x0,d)
            % solve first fail, so solve the second one to make sure the
            % continuous u F
            
            if d < 3 % transDSP transSSP
                obj.H = blkdiag(zeros(obj.n_u),eye(obj.n_rod),eye(obj.n_GRF),zeros(obj.hol_size_2D), 0);
                obj.g = -[zeros(obj.n_u,1);2*ones(obj.n_rod + obj.n_GRF,1);0].*x0;
            else
                obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_hol),0.*eye(obj.n_GRF),zeros(obj.hol_size_2D),obj.p);
                obj.g = 0.*obj.g + zeros(obj.n_u_hat+1,1);
            end
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 4000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g,A,obj.lb,obj.ub,lbA,ubA,options);

%             options = optimoptions('quadprog','Display','off');
%             [xsol,fval,exitflag, numiter] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,x0,options);
%             
            u = xsol(1:obj.n_u);
            Frod = xsol(obj.n_u+1:obj.n_u+obj.n_rod);
            F_GRF = xsol(obj.n_u+obj.n_rod+1:obj.n_u+obj.n_rod+obj.n_GRF);
            delta = xsol(end);
        end
        
        function [u,exitflag,Frod,F_GRF,delta, fval, numiter] = safesolve_withFhol(obj,x0,d)
            obj.RemoveConstraints_withFhol(); 
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 4000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g,A,obj.lb,obj.ub,lbA,ubA,options);
            
            if exitflag ~= 0
                options = optimoptions('quadprog','Display','off');
                [xsol,fval,exitflag, output] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,x0,options);
                numiter = output.iterations;
                exitflag = obj.convertMatlabQuadprog2QPOASES(exitflag);
            end
            
            u = xsol(1:obj.n_u);
            Frod = xsol(obj.n_u+1:obj.n_u+obj.n_rod);
            F_GRF = xsol(obj.n_u+obj.n_rod+1:obj.n_u+obj.n_rod+obj.n_GRF);
            delta = xsol(end);
        end
       
    end
    
    methods % Operational Space for robots
              
        function ConfigMidFootQP_OSC(obj)
            obj.lb = -obj.motorTorqueLimits';
            obj.ub = -obj.lb;
%             
            obj.n_u = 4; 
            obj.n_GRF = 2;

            obj.Aeq = [];
            obj.Beq = [];
            obj.Aneq = [];
            obj.Bneq = [];
            % friction cone u_friction_cone_RightMidFoot(RightForceXYZ, miu)>=0;
            % moment constraint; not consider yaw moment for now.
            appx_miu = obj.miu; % 
            
            obj.ContactMatrix = [1 -appx_miu;
                                -1 -appx_miu;
                                0 -1];
            
            rowsGRF = size(obj.ContactMatrix,1);
            
            % GRF force Fznormal>0 frictionxy<uFnorm, My <lt lh*Fz Mz <lt lh * Fy
            obj.AneqGRF = obj.ContactMatrix;
            obj.BneqGRF = zeros(rowsGRF,1);
            
        end
        
        function configureCost_OSC(obj, outputs, eom)
            %%% min||Au - b||
            J = outputs.RD2.J_y;
            dJ = outputs.RD2.dJ_y;
            y = outputs.RD2.y;
            dy = outputs.RD2.dy;
            
         %   KP_SSP =  blkdiag(8000, 5000, 15000, 5000);
         %   KD_SSP =  blkdiag(800, 500, 1500, 500);
            
            KP_SSP =  blkdiag(1500, 500, 1500, 500);
            KD_SSP =  blkdiag(150, 50, 150, 50);
            ddy_des = - KP_SSP*y - KD_SSP*dy;
            
            A = J*eom.gx_hat_bottom;
            b = J*eom.fx_hat_bottom + dJ*eom.dq - ddy_des;
            %%% J*fx too large than ddy_des
            W = 1; %1e+4;
            L = 1e-5; % prevent singular
            obj.H = W*A'*A + L*eye(4);
            obj.g = (W*b'*A)';
        end
        
        function ConfigureConstraints_OSC(obj, eom, Contact)
            obj.contact  = Contact;
            A_GRF = [obj.AneqGRF*eom.A_Fhol2U];
            b_GRF = -obj.AneqGRF*eom.b_Fhol2U;
            obj.Aneq = A_GRF;
            obj.Bneq = b_GRF;
        end
       
        function RemoveConstraints_OSC(obj, eom)
            %%% 
            disp('something wrong:: want to fly?')
        end
       
        function [u, exitflag, F_GRF, fval, numiter] = solve_OSC(obj, eom)
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 10000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            [xsol,fval,exitflag, numiter] = qpOASES(obj.H, obj.g, A, obj.lb, obj.ub, lbA, ubA, options);
%             options = optimoptions('quadprog','Display','off');
%             [xsol,fval,exitflag, numiter] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,[],options);
            u = xsol;
            F_GRF = eom.A_Fhol2U*u + eom.b_Fhol2U; % if flight phase, its rod force
        end
        
        function [u, exitflag, F_GRF, fval, numiter] = safesolve_OSC(obj, x0, eom)
            obj.RemoveConstraints_OSC(eom); 
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 4000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            % [xsol,fval,exitflag, numiter] = qpOASES(obj.H, obj.g, A, [-Inf;-Inf;-Inf;-Inf],  [Inf;Inf;Inf;Inf], lbA, ubA, options);

            [xsol,fval,exitflag, numiter] = qpOASES(obj.H, obj.g, A, obj.lb, obj.ub, lbA, ubA, options);
            
            if exitflag ~= 0
                options = optimoptions('quadprog','Display','off');
                [xsol,fval,exitflag, output] = quadprog(obj.H, obj.g, obj.Aneq, obj.Bneq, obj.Aeq, obj.Beq, obj.lb, obj.ub, x0, options);
                numiter = output.iterations; 
                exitflag = obj.convertMatlabQuadprog2QPOASES(exitflag);
            end
            if isempty(xsol) 
                xsol = zeros(4,1); 
                disp('quadprog failed')
            end
            u = xsol;
            F_GRF = eom.A_Fhol2U*u + eom.b_Fhol2U;
        end
        
     end
    
    methods % QP that only uses torques and relaxation variable as optimization variables
        
        function ConfigMidFootQP_onlyU(obj)
            obj.lb = [-obj.motorTorqueLimits;-inf];
            obj.ub = -obj.lb;
            obj.n_u = 10; obj.n_rod = 2; 
            obj.n_GRF = 10;
            obj.n_delta = 1; 
            if strcmp(obj.model2D3D, '2D')
                obj.hol_size_2D = 3;
            else
                obj.hol_size_2D = 0;
            end
            n_arguments = [obj.n_u, obj.n_delta];
            obj.n_total = sum(n_arguments); 
            obj.n_u_hat = obj.n_total - 1; 
            obj.Aeq = [];
            obj.Beq = [];
            obj.Aneq = [];
            obj.Bneq = [];
            % friction cone u_friction_cone_RightMidFoot(RightForceXYZ, miu)>=0;
            % moment constraint; not consider yaw moment for now.
            %appx_miu = sqrt(2)/2*obj.miu; % polyhyd appx
            appx_miu = obj.miu; % 
            
            if strcmp(obj.model2D3D ,'3D')
                obj.ContactMatrix = [
                    1 0 -appx_miu 0 0;
                    -1 0 -appx_miu 0 0;   % forward friction
                    0 1 -appx_miu 0 0;
                    0 -1 -appx_miu 0 0;   % lateral friction
                    0 0 -1 0 0;           % normal force
                    0 0 -obj.half_foot_length 1 0;
                    0 0 -obj.half_foot_length -1 0;  %forward zmp
                    0 0 -obj.half_foot_length*obj.miu 0 1;
                    0 0 -obj.half_foot_length*obj.miu 0 -1;
                    ];
            elseif strcmp(obj.model2D3D ,'2D')
                obj.ContactMatrix = [
                 %   1 0 -appx_miu 0 0; % due to the holonomic constraint   on the pelvis, this friction has to be commented out
                 %  -1 0 -appx_miu 0 0;   % forward friction % this needs to be removed due to  
                    0 0 -1 0 0;           % normal force
                 0 0 -obj.half_foot_length 1 0;
                   0 0 -obj.half_foot_length -1 0;  %forward zmp
                    ];
            elseif strcmp(obj.model2D3D ,'2Dsoft')
                obj.ContactMatrix = [
             %       1 0 -appx_miu 0 0;
               %    -1 0 -appx_miu 0 0;   % forward friction
                    0 0 -1 0 0;           % normal force
                    0 0 -obj.half_foot_length 1 0;
                    0 0 -obj.half_foot_length -1 0;  %forward zmp
                    ];
            end
            rowsGRF = size(obj.ContactMatrix,1);
            
            % GRF force Fznormal>0 frictionxy<uFnorm, My <lt lh*Fz Mz <lt lh * Fy
            obj.AneqGRF = [zeros(rowsGRF,obj.n_rod + obj.n_GRF/2), obj.ContactMatrix, zeros(rowsGRF, obj.hol_size_2D);
                           zeros(rowsGRF,obj.n_rod), obj.ContactMatrix, zeros(rowsGRF,obj.n_GRF/2 + obj.hol_size_2D)];
            obj.BneqGRF = zeros(2*rowsGRF,1);
            
            obj.AneqGRFleft = [zeros(rowsGRF, obj.n_rod),obj.ContactMatrix,zeros(rowsGRF, obj.hol_size_2D)];
            obj.BneqGRFleft = zeros(rowsGRF,1);
            obj.AneqGRFright = [zeros(rowsGRF, obj.n_rod),obj.ContactMatrix,zeros(rowsGRF,obj.hol_size_2D)];
            obj.BneqGRFright = zeros(rowsGRF,1);
        end
        
        function configureCost_onlyU(obj, A_hat, QP_Lf, smoothingQP, contact,d, eom)
            ATA = A_hat'*A_hat;
            obj.H = blkdiag(ATA, obj.p);
            obj.g = [(QP_Lf'*A_hat)';zeros(obj.n_delta, 1)];
            
            grfMin = blkdiag(1,1,0,1,1);
            grfMIN = blkdiag(grfMin, grfMin);
            obj.H =  obj.H + smoothingQP*blkdiag(eye(obj.n_u), 0); %smoother and faster
            %obj.H =  blkdiag(eye(obj.n_u),zeros(obj.n_hol),eye(obj.n_GRF),obj.p);
            % if sum(desiredGRF)~= 0
            %obj.H = obj.H + blkdiag(zeros(obj.n_u),obj.P_GRF,zeros(obj.n_hol+1));
            % obj.g = obj.g - [zeros(obj.n_u,1);obj.P_GRF*desiredGRF;zeros(obj.n_hol+1,1)];
            % end
%             obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_hol),0.*eye(obj.n_GRF),obj.p);
%             obj.g = 0.*obj.g + zeros(obj.n_u_hat+1,1);
            
            % specific for walking
            % if minimize u; this could work
            if ( strcmp(contact,'DoubleSupportFromLeft')||strcmp(contact,'DoubleSupportFromRight')) 
                if d == 1
                    obj.H = 1.*obj.H + 0*blkdiag(eye(10), obj.p);
                    obj.g = 1.*obj.g + zeros(obj.n_u_hat+1,1);
%                     obj.H = blkdiag(eom.A_Fhol2U'*eom.A_Fhol2U, obj.p); 
%                     obj.g = [2*(eom.b_Fhol2U'*eom.A_Fhol2U)';zeros(1, 1)];
                else
                    obj.H = .0*obj.H + 1*blkdiag(eye(10), obj.p);
                    obj.g = .0*obj.g + zeros(obj.n_u_hat+1,1);
                    % minimize Fhol 
                   % obj.H = blkdiag(eom.A_Fhol2U'*eom.A_Fhol2U, obj.p); 
                   % obj.g = [2*(eom.b_Fhol2U'*eom.A_Fhol2U)';zeros(1, 1)];
                end
            end
%             obj.H = .0*obj.H + 1*blkdiag(eye(10),obj.p);
%             obj.g = .0*obj.g + zeros(obj.n_u_hat+1,1);
        end
        
        function ConfigureConstraints_onlyU(obj, A_hat,QP_Lf,psi0,psi1,eom, Contact)
            obj.contact  = Contact; 
            %% Inequality
            % CLF: psi0 + psi1'*(A_hat*u_hat + Lf) <= delta
            obj.Aclf = [psi1'*A_hat, -1];
            obj.Bclf = -psi1'*QP_Lf - psi0;
            
            switch Contact
                case AmberConstants.RightFootContact
                    A_GRF = [obj.AneqGRFright*eom.A_Fhol2U, zeros(size(obj.AneqGRFright,1),1)];
                    b_GRF = -obj.AneqGRFright*eom.b_Fhol2U;
                case AmberConstants.LeftFootContact
                    A_GRF = [obj.AneqGRFleft*eom.A_Fhol2U, zeros(size(obj.AneqGRFleft,1),1)];
                    b_GRF = -obj.AneqGRFleft*eom.b_Fhol2U;
                case 'Flight'
                    A_GRF = [];
                    b_GRF = [];
                otherwise %  'DoubleSupportFromRight' 'DoubleSupportFromLeft' 'Prejump' 'landing' or something% DSP is the same
                    A_GRF = [obj.AneqGRF*eom.A_Fhol2U, zeros(size(obj.AneqGRF,1),1)];
                    b_GRF = -obj.AneqGRF*eom.b_Fhol2U;   
            end
            obj.Aneq = [obj.Aclf; A_GRF];
            obj.Bneq = [obj.Bclf; b_GRF];
            obj.Aeq = [];
            obj.Beq = [];
        end
       
        function RemoveConstraints_onlyU(obj, eom)
            switch obj.contact
                case 'DoubleSupportFromRight' % next domain is left foot contact
                    xinv = inv(eom.J_hol_leftContact*(eom.Mass\eom.J_hol_leftContact'));
                    % Fhol = A_Fhol2U*u + b_Fhol2U
                    A_Fhol2U = -xinv*eom.J_hol_leftContact*(eom.Mass\eom.Bu);
                    b_Fhol2U = xinv*(eom.J_hol_leftContact*(eom.Mass\eom.h) - eom.dJ_hol_leftContact*eom.dq);
                    
                    A_GRF = [obj.AneqGRFleft*A_Fhol2U,  zeros(size(obj.AneqGRFleft,1),1)];
                    b_GRF = -obj.AneqGRFleft*b_Fhol2U;
                    obj.Aneq = [obj.Aclf;A_GRF];
                    obj.Bneq = [obj.Bclf;b_GRF];
                case 'DoubleSupportFromLeft'
                    xinv = inv(eom.J_hol_rightContact*(eom.Mass\eom.J_hol_rightContact'));
                    % Fhol = A_Fhol2U*u + b_Fhol2U
                    A_Fhol2U = -xinv*eom.J_hol_rightContact*(eom.Mass\eom.Bu);
                    b_Fhol2U = xinv*(eom.J_hol_rightContact*(eom.Mass\eom.h) - eom.dJ_hol_rightContact*eom.dq);
                    
                    A_GRF = [obj.AneqGRFright*A_Fhol2U,  zeros(size(obj.AneqGRFright,1),1)];
                    b_GRF = -obj.AneqGRFright*b_Fhol2U;
                    obj.Aneq = [obj.Aclf;A_GRF];
                    obj.Bneq = [obj.Bclf;b_GRF];
                case 'PreJump'
                    % holonomic constraints on rod only, no force limit
                    % TO CONSIDER, the accelereation of the foot z is only
                    % admissible to be up/positive, which is a constraint
                    % on u, JD-1B u - JD-1h + Jdot qdot >0
                    obj.Aneq = obj.Aclf;
                    obj.Bneq = obj.Bclf;
                case 'PeriodicStance' % the same to Prejump
                    obj.Aneq = obj.Aclf;
                    obj.Bneq = obj.Bclf;
                otherwise
                    disp('QP fails in SSP')
            end
            
        end
        
        function [u,exitflag,Frod,F_GRF,delta, fval, numiter] = solve_onlyU(obj, eom)
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 10000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g, A,obj.lb,obj.ub,lbA,ubA,options);
%             options = optimoptions('quadprog','Display','off');
%             [xsol,fval,exitflag, numiter] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,[],options);
            u = xsol(1:obj.n_u);
            F_hol = eom.A_Fhol2U*u + eom.b_Fhol2U; % if flight phase, its rod force
            
            Frod = F_hol(1:2); 
            F_GRF = zeros(10,1); %F_hol(3:end-3); 
            delta = xsol(end);
        end
        
        function [u,exitflag,Frod,F_GRF,delta, fval, numiter] = safesolve_onlyU(obj,x0,d, eom)
            obj.RemoveConstraints_onlyU(eom); 
            [A, lbA, ubA] = obj.convertConstraints();
            nWSR = 4000;
            options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g,A,obj.lb,obj.ub,lbA,ubA,options);
            
            if exitflag ~= 0
                options = optimoptions('quadprog','Display','off');
                [xsol,fval,exitflag, output] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,x0,options);
                numiter = output.iterations;
                exitflag = obj.convertMatlabQuadprog2QPOASES(exitflag);
            end
            if isempty(xsol) 
                xsol = zeros(11,1); 
                disp('quadprog failed')
            end
            u = xsol(1:obj.n_u);
            F_hol = eom.A_Fhol2U*u + eom.b_Fhol2U;
            
            Frod = F_hol(1:2);
            F_GRF = zeros(10,1);%F_hol(3:12);
            delta = xsol(end);
        end
        
    end
 
    methods % QP that also takes ddq into variable,: QPoasas has difficulty to find solution, works with quadprog 
        function ConfigMidFootQP_withddq(obj)
            obj.n_u = 10; obj.n_rod = 2; 
            obj.n_GRF = 10;
            obj.n_delta = 1; 
            obj.ddq_size = 22; 
            if strcmp(obj.model2D3D, '2D')
                obj.hol_size_2D = 3;
            else
                obj.hol_size_2D = 0;
            end
            n_arguments = [obj.n_u, obj.n_rod, obj.n_GRF, obj.hol_size_2D, obj.ddq_size, obj.n_delta];
            obj.n_total = sum(n_arguments); 
            
            obj.n_u_hat = obj.n_total - obj.ddq_size - obj.n_delta; 
            
            obj.lb = [-obj.motorTorqueLimits;-inf*ones(obj.n_total-obj.n_u,1)];
            obj.ub = -obj.lb;

            appx_miu = obj.miu; % 
            if strcmp(obj.model2D3D ,'3D')
                obj.ContactMatrix = [
                    1 0 -appx_miu 0 0;
                    -1 0 -appx_miu 0 0;   % forward friction
                    0 1 -appx_miu 0 0;
                    0 -1 -appx_miu 0 0;   % lateral friction
                    0 0 -1 0 0;           % normal force
                    0 0 -obj.half_foot_length 1 0;
                    0 0 -obj.half_foot_length -1 0;  %forward zmp
                    0 0 -obj.half_foot_length*obj.miu 0 1;
                    0 0 -obj.half_foot_length*obj.miu 0 -1;
                    ];
            elseif strcmp(obj.model2D3D ,'2D')
                obj.ContactMatrix = [
                 %   1 0 -appx_miu 0 0; % due to the holonomic constraint   on the pelvis, this friction has to be commented out
                 %  -1 0 -appx_miu 0 0;   % forward friction % this needs to be removed due to  
                    0 0 -1 0 0;           % normal force
                  0 0 -obj.half_foot_length 1 0;
                   0 0 -obj.half_foot_length -1 0;  %forward zmp
                    ];
            elseif strcmp(obj.model2D3D ,'2Dsoft')
                obj.ContactMatrix = [
             %       1 0 -appx_miu 0 0;
               %    -1 0 -appx_miu 0 0;   % forward friction
                    0 0 -1 0 0;           % normal force
                    0 0 -obj.half_foot_length 1 0;
                    0 0 -obj.half_foot_length -1 0;  %forward zmp
                    ];
            end
            rowsGRF = size(obj.ContactMatrix,1);

            % GRF force Fznormal>0 frictionxy<uFnorm, My <lt lh*Fz Mz <lt lh * Fy
            obj.AneqGRF = [zeros(rowsGRF,obj.n_u + obj.n_rod), obj.ContactMatrix, zeros(rowsGRF,obj.n_GRF/2 + obj.hol_size_2D + obj.ddq_size + obj.n_delta);
                zeros(rowsGRF,obj.n_u + obj.n_rod + obj.n_GRF/2), obj.ContactMatrix, zeros(rowsGRF,obj.hol_size_2D + obj.ddq_size +obj.n_delta)];
            obj.BneqGRF = zeros(2*rowsGRF,1);
            
            obj.AneqGRFleft = [zeros(rowsGRF,obj.n_u + obj.n_rod), obj.ContactMatrix, zeros(rowsGRF,obj.n_GRF/2 + obj.hol_size_2D + obj.ddq_size + obj.n_delta)];
            obj.BneqGRFleft = zeros(rowsGRF,1);
            obj.AneqGRFright = [zeros(rowsGRF,obj.n_u + obj.n_rod + obj.n_GRF/2), obj.ContactMatrix, zeros(rowsGRF,obj.hol_size_2D + obj.ddq_size +obj.n_delta)];
            obj.BneqGRFright = zeros(rowsGRF,1);
            
            obj.AeqGRFleft =  [zeros(obj.n_GRF/2, obj.n_u + obj.n_rod),eye(obj.n_GRF/2),zeros(obj.n_GRF/2,obj.n_GRF/2+ obj.hol_size_2D + obj.ddq_size + obj.n_delta)];
            obj.BeqGRFleft = zeros(obj.n_GRF/2,1);
            obj.AeqGRFright = [zeros(obj.n_GRF/2, obj.n_u + obj.n_rod + obj.n_GRF/2),eye(obj.n_GRF/2),zeros(obj.n_GRF/2,obj.hol_size_2D + obj.ddq_size +obj.n_delta)];
            obj.BeqGRFright = zeros(obj.n_GRF/2,1);
        end
        
        function configureCost_withddq(obj, A_hat, QP_Lf, smoothingQP,contact,d)
            ATA = A_hat'*A_hat; 
            ddq_weights = 0.; 
            obj.H = blkdiag( ATA, ddq_weights*ones(obj.ddq_size) , obj.p );
            obj.g = [ 2*(QP_Lf'*A_hat)'; zeros(obj.ddq_size + obj.n_delta, 1)];
            
            grfMin = blkdiag(1,1,0,1,1);
            grfMIN = blkdiag(grfMin, grfMin);
            obj.H =  obj.H + smoothingQP*blkdiag(ones(obj.n_u),zeros(obj.n_rod),0*grfMIN,zeros(obj.ddq_size + obj.hol_size_2D), 0); %smoother and faster
            % if minimize u; this could work
            if ( strcmp(contact,'DoubleSupportFromLeft')||strcmp(contact,'DoubleSupportFromRight'))  && d> 1
                obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_rod),0.*eye(obj.n_GRF + obj.hol_size_2D), zeros(obj.ddq_size), obj.p);
                obj.g = 0.*obj.g + zeros(obj.n_u_hat + obj.ddq_size +1,1);
            end
        end
        
        function ConfigureConstraints_withddq(obj, A_hat, QP_Lf, psi0, psi1, eom, Contact)
            obj.contact = Contact; 
            %% equality constraints
            % B_hat*u_hat + D*ddq +h = 0;
            obj.Beq_eom = eom.CGvec - eom.springVec;
            obj.Aeq_eom = [eom.B_hat, -eom.Mass, zeros(length(obj.Beq_eom),obj.n_delta)];
            
            % holonomic rod
            obj.Beq_rod = - eom.dJ_rod*eom.dq;
            obj.Aeq_rod = [zeros(length(obj.Beq_rod),obj.n_u + obj.n_rod + obj.n_GRF + obj.hol_size_2D), eom.J_rod, zeros(length(obj.Beq_rod),obj.n_delta)];
            
            % holonomic 2D
            if strcmp(obj.model2D3D, '2D')
                obj.Beq_2D = - eom.dJ_2Dpelvis*eom.dq;
                obj.Aeq_2D = [zeros(length(obj.Beq_2D),obj.n_u + obj.n_rod + obj.n_GRF + obj.hol_size_2D), eom.J_2Dpelvis, zeros(length(obj.Beq_2D),obj.n_delta)];
            else
                obj.Beq_2D = [];
                obj.Aeq_2D = [];
            end
            % left contact hol 
            obj.Beq_left = - eom.dJ_leftStance*eom.dq;
            obj.Aeq_left = [zeros(length(obj.Beq_left),obj.n_u + obj.n_rod + obj.n_GRF + obj.hol_size_2D), eom.J_leftStance, zeros(length(obj.Beq_left),obj.n_delta)];
            
             % right contact hol 
            obj.Beq_right = - eom.dJ_rightStance*eom.dq;
            obj.Aeq_right = [zeros(length(obj.Beq_right),obj.n_u + obj.n_rod + obj.n_GRF + obj.hol_size_2D), eom.J_rightStance, zeros(length(obj.Beq_right),obj.n_delta)];
            
            %% Inequality
            % CLF: psi0 + psi1'*(A_hat*u_hat + Lf) <= delta
            obj.Aclf = [psi1'*A_hat, zeros(1, obj.ddq_size), -1];
            obj.Bclf = -psi1'*QP_Lf - psi0;
            
            switch Contact
                case AmberConstants.RightFootContact
                    obj.Aneq = [obj.Aclf; obj.AneqGRFright];
                    obj.Bneq = [obj.Bclf; obj.BneqGRFright];
                    obj.Aeq = [obj.Aeq_eom;  obj.AeqGRFleft; obj.Aeq_rod; obj.Aeq_2D; obj.Aeq_right];
                    obj.Beq = [obj.Beq_eom;  obj.BeqGRFleft; obj.Beq_rod; obj.Beq_2D; obj.Beq_right];
                case AmberConstants.LeftFootContact
                    obj.Aneq = [obj.Aclf; obj.AneqGRFleft];
                    obj.Bneq = [obj.Bclf; obj.BneqGRFleft];
                    obj.Aeq = [obj.Aeq_eom;  obj.AeqGRFright; obj.Aeq_rod; obj.Aeq_2D; obj.Aeq_left];
                    obj.Beq = [obj.Beq_eom;  obj.BeqGRFright; obj.Beq_rod; obj.Beq_2D; obj.Beq_left];
                case 'DoubleSupportFromRight' % DSP is the same
                    obj.Aneq = [obj.Aclf; obj.AneqGRF];
                    obj.Bneq = [obj.Bclf; obj.BneqGRF];
                    obj.Aeq = [obj.Aeq_eom; obj.Aeq_rod; obj.Aeq_2D; obj.Aeq_left; obj.Aeq_right];
                    obj.Beq = [obj.Beq_eom; obj.Beq_rod; obj.Beq_2D; obj.Beq_left; obj.Beq_right];
                case 'DoubleSupportFromLeft'
                    obj.Aneq = [obj.Aclf; obj.AneqGRF];
                    obj.Bneq = [obj.Bclf; obj.BneqGRF];
                    obj.Aeq = [obj.Aeq_eom; obj.Aeq_rod; obj.Aeq_2D; obj.Aeq_left; obj.Aeq_right];
                    obj.Beq = [obj.Beq_eom; obj.Beq_rod; obj.Beq_2D; obj.Beq_left; obj.Beq_right];
            end
        end
        
        function [A, lbA, ubA] = convertConstraints_withddq(obj)
            A = [obj.Aneq; obj.Aeq];
            ubA = [obj.Bneq; obj.Beq];
            lbA = [-inf*ones(size(obj.Bneq));obj.Beq];
        end
                
        function RemoveConstraints_withddq(obj)
            switch obj.contact
                case 'DoubleSupportFromRight' % the next domain is left foot contact
                    obj.Aneq = [obj.Aclf; obj.AneqGRFleft];
                    obj.Bneq = [obj.Bclf; obj.BneqGRFleft];
                    obj.Aeq = [obj.Aeq_eom;  obj.AeqGRFright; obj.Aeq_rod; obj.Aeq_2D; obj.Aeq_left];
                    obj.Beq = [obj.Beq_eom;  obj.BeqGRFright; obj.Beq_rod; obj.Beq_2D; obj.Beq_left];
                case 'DoubleSupportFromLeft' % the next domain is right foot contact
                     obj.Aneq = [obj.Aclf; obj.AneqGRFright];
                    obj.Bneq = [obj.Bclf; obj.BneqGRFright];
                    obj.Aeq = [obj.Aeq_eom;  obj.AeqGRFleft; obj.Aeq_rod; obj.Aeq_2D; obj.Aeq_right];
                    obj.Beq = [obj.Beq_eom;  obj.BeqGRFleft; obj.Beq_rod; obj.Beq_2D; obj.Beq_right];
                otherwise
                    disp('QP fails in SSP');
            end
        end
           
        function [u, exitflag, Frod, F_GRF, ddq, delta, fval, numiter] = solve_withddq(obj)
            %             [A, lbA, ubA] = obj.convertConstraintsNew();
            %             nWSR = 10000;
            %             options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            %             [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g,A,obj.lb,obj.ub,lbA,ubA,options);
            %
            options = optimoptions('quadprog','Display','off');
            [xsol,fval,exitflag, output] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,obj.initialGuess,options);
            if exitflag ~= 1
                disp('QP fail');
            end
            % exitflag = 1; converged; 
            % 0 : number of iterations exceeded; 
            % -2 : infeasible 
            numiter = output.iterations;
            u = xsol(1:obj.n_u);
            Frod = xsol(obj.n_u+1:obj.n_u+obj.n_rod);
            F_GRF = xsol(obj.n_u+obj.n_rod+1:obj.n_u+obj.n_rod+obj.n_GRF);
            ddq = xsol( (end-obj.ddq_size):(end-1));
            delta = xsol(end);
            obj.initialGuess = xsol; 
        end
        
        function [u, exitflag, Frod, F_GRF, ddq, delta, fval, numiter] = safesolve_withddq(obj,x0,d)
            % solve first fail, so solve the second one to make sure the
            % continuous u F
            %             if d < 3 % transDSP transSSP
            %                 obj.H = blkdiag(zeros(obj.n_u),eye(obj.n_rod),eye(obj.n_GRF),eye(obj.hol_size_2D), eye(obj.ddq_size), 0);
            %                 obj.g = -[zeros(obj.n_u,1);2*ones(obj.n_rod + obj.n_GRF,1);0].*x0;
            %             else
            %% QP fails at the end of DSP
            %             obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_rod),0.*eye(obj.n_GRF + obj.hol_size_2D + obj.ddq_size),obj.p);
            %             obj.g = 0.*obj.g + zeros(obj.n_u_hat+ obj.ddq_size +1,1);
            %             options = qpOASES_options('maxIter',nWSR,'printLevel',1);
            %             [xsol,fval,exitflag, numiter] = qpOASES(obj.H,obj.g,A,obj.lb,obj.ub,lbA,ubA,options);
                        obj.H = 0.*obj.H + 1*blkdiag(eye(10),0*eye(obj.n_rod),0.*eye(obj.n_GRF + obj.hol_size_2D),1*eye(obj.ddq_size),obj.p);
                        obj.g = 0.*obj.g + zeros(obj.n_u_hat+ obj.ddq_size +1,1);
            
            % obj.RemoveConstraintsNew();
            options = optimoptions('quadprog','Display','off');
            [xsol,fval,exitflag, output] = quadprog(obj.H,obj.g,obj.Aneq,obj.Bneq,obj.Aeq,obj.Beq,obj.lb,obj.ub,x0,options);
            if exitflag ~= 1
                disp('QP fail');
            end
            numiter = output.iterations;
            u = xsol(1:obj.n_u);
            Frod = xsol(obj.n_u+1:obj.n_u+obj.n_rod);
            F_GRF = xsol(obj.n_u+obj.n_rod+1:obj.n_u+obj.n_rod+obj.n_GRF);
            ddq = xsol( (end-obj.ddq_size):(end-1));
            
            delta = xsol(end);
        end
    end
        
end

