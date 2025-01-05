classdef GenOutputsFiveLinkWalker < handle
    
    properties
        % all possible outputs to choose
        % rd = 2  y = a(q); y_dot = J*dq; y_ddot = J*ddq + dJ*dq;
        pelvisPos
        comPos
        comVel
        
        RD2
        RD2unscaled
        RD1
        
        legSelector
        % time based save
        q
        dq
        
        % phase
        phaseSave   = struct('t',[], 'value', []); % for SSP
        LIP_DSP_dy0 = 0;
        motorIdx = 4:7;
        
        %%%% virtual constraints %%%
        a
        p
        p0
        pN
        polyMat
        
        Jmotors = zeros(4, 7);
        dJmotors = zeros(4, 7);
        tvec
        qmdes
        dqmdes
        timeOrPhase = 'phase'
        
        gaitType %%%% 'P1', 'P2'
        %%%%% walking construction via LIP %%%%%%%%%
        impactVel %%%  impact velocity 0.05
        swingHeight %%% swing foot height
        swingZtraj
        tSSP
        massHeight
        g = 9.81;
        Ahat0
        Bhat0
        vxdes
        Xdes
        XLdes
        XRdes
        KLQR
        Kdeadbeat
        KposDeadbeat
        period
        smfVec
        smfVecVel
        smfVec2
        smfVecVel2
        SSPnsf0 = [];
        SSPnsfVel0 = [];
        stanceToePos
        qPelvis0 = 0; 
        Nmovingavg = 50;
        
        %%% LIP stuff
        accCOMvec
        velCOMvec
        pVec = []
        vVec = []
        uVec = []; 
        stepL_left = 0 %%% one step size for P2 orbits
        d2
        
        %%% DM control 
        y_kp1 = zeros(4,1); 
    end
    
    methods
        function obj = GenOutputsFiveLinkWalker(logOpt, tSSP, COMheight, swingHeight, impactVel, vxdes, gaitType,stepL_left, period)
            obj.tSSP = tSSP;
            obj.massHeight = COMheight ;
            obj.swingHeight = swingHeight;
            obj.impactVel = impactVel ;
            obj.vxdes = vxdes;
            obj.period = period;
            obj.gaitType = gaitType;
            obj.stepL_left = stepL_left;
            a = logOpt.alpha;
            p = logOpt.p;
            obj.RD2 = [];
            obj.RD1 = struct('y',[],'Aq',[],'JA',[]);
            obj.legSelector = eye(2);
            
            for i = 4:7
                obj.Jmotors(i-3, i) = 1;
            end
            
            obj.a = reshape(a,4, []);
            obj.p = p;
            % obj.VRreconstruct();
            % obj.timeBasedConstruct(logOpt);
            % obj.VRangleConstruct(logOpt);
            obj.SwingFootZconstruct();
            obj.initializeLIP();
            
            obj.y_kp1 = [ 0; 0;
                obj.massHeight;
                obj.qPelvis0];
            obj.accCOMvec = zeros(1, obj.Nmovingavg);
            obj.velCOMvec = zeros(1, obj.Nmovingavg);
        end
        %%%
        
        function initializeLIP(obj)
            TS = obj.tSSP;
            TD = 0;
            lambda = sqrt(obj.g/obj.massHeight);
            eATs = [ cosh( TS*lambda), sinh(TS*lambda)/lambda;
                lambda*sinh(TS*lambda), cosh(TS*lambda)];
            
            %%%
            obj.Bhat0 = eATs*[-1;
                0];
            obj.Ahat0 = eATs*[1,  TD;
                0, 1];
            
            sigma_En = lambda*tanh(TS/2*lambda);
            sigma_Ep = lambda*coth(TS/2*lambda);
            
            obj.Xdes = [obj.vxdes/sigma_Ep; obj.vxdes];
            %%%%% for P2 orbits %%%%%%%%%
            % characteristic line dotX = k x + d;
            d = lambda^2*(sech(lambda*TS/2))^2*obj.vxdes*(TS + TD)/(lambda^2*TD + 2*sigma_En);
            % select the boundary condition. % x1_F
            uL = obj.stepL_left;
            TD = 0;
            pLdes = (uL - TD*d)/(2 + TD*sigma_En); % m  custom selection
            vLdes = sigma_En*pLdes + d;
            [pRdes, vRdes, ~] = obj.LIP_sol(lambda, -pLdes, vLdes,TS);
            obj.XLdes = [pLdes; vLdes];
            obj.XRdes = [pRdes; vRdes];
            
            %%% update dLQR gain
            %%%% LQR %%%%%%%%%%%%%%%%%%
            N = zeros(2,1);
            R = 1;
            Q = 10*eye(2);
            [obj.KLQR, ~, ~] = dlqr(obj.Ahat0, obj.Bhat0, Q, R, N);
            obj.KLQR = -obj.KLQR; %%% sign convention
            %%%% deadbeat gain %%%%%%% 
            obj.Kdeadbeat = [1, TD + coth(TS*lambda)/lambda] ;
            obj.KposDeadbeat =  [1, TD + tanh(TS*lambda)/lambda] ;
            
            swingConstructT = 0:obj.period:3;
            obj.smfVec = smf(swingConstructT, [obj.tSSP/10 obj.tSSP*9/10]);
            obj.smfVecVel = [diff(obj.smfVec)/obj.period, 0];
            obj.smfVec2 = smf(swingConstructT, [0 obj.tSSP/5]);
            obj.smfVecVel2 = [diff(obj.smfVec2)/obj.period, 0];
        end
    end
    
    %%%%%% note the selection of outputs matters of-course
    %%%%%% plevis angle; leg length; swing foot z;
    methods %(SSP)
        function obj = getOutputs(obj, q, dq, ddq, tIdx, Contact)
            OutputScaling = struct('comXYZ', [1; 1; 1], 'pelvisPitch', 1, 'SwingToeXZ', [10;1]);
            obj.SSP(q, dq, ddq, tIdx, OutputScaling, Contact);
        end

        % use real COM, encoded by VirtualLeg Angle
        function obj = SSP(obj, q, dq, ddq, tIdx, OutputScaling, Contact)
            %  COMz, swingZ, swingX, pelvis
            if isempty(obj.SSPnsf0)
                if  Contact == AmberConstants.LeftFootContact
                    obj.SSPnsf0      = pRightToe(q);
                    obj.SSPnsfVel0   = vRightToe([q;dq]);
                    obj.stanceToePos = pLeftToe(q);
                elseif  Contact == AmberConstants.RightFootContact
                    obj.SSPnsf0      = pLeftToe(q);
                    obj.SSPnsfVel0   = vLeftToe([q;dq]);
                    obj.stanceToePos = pRightToe(q);
                end
            end
            obj.getCOM(q, dq, ddq);
            % calculate outputs Pos vel J dJ
            SwingFootPos  = obj.OutputSwingFootPos(q, dq, OutputScaling, Contact);
            yPelvisEuler  = obj.OutputPelvisEuler(q, dq, OutputScaling) ;
            yCOMpos       = obj.OutputCOMpos(q, dq, OutputScaling) ;
            
            % swing foot planning
            [swingX, stepLengthDes, stepLengthVelDes] = obj.swingFootPlanning(q,dq, tIdx, Contact);
            
            obj.comPos = yCOMpos.PosOri(1);
            obj.comVel = yCOMpos.VelOri(1);
            N  = length(obj.swingZtraj.pos);
            tIdx(tIdx>N) = N;
            swingZ = struct('pos', obj.swingZtraj.pos(tIdx), 'vel', obj.swingZtraj.vel(tIdx));
            %   swingZ =  struct('pos', 0, 'vel', 0);
            %   swingX =  struct('pos', 0, 'vel', 0);
            %
            comZ =    struct('pos', obj.massHeight, 'vel', 0);
            qPelvis = struct('pos', obj.qPelvis0,   'vel', 0);
            % desired outputs
            y_des_unscaled =  [ swingX.pos+obj.stanceToePos(1); swingZ.pos;     comZ.pos;  qPelvis.pos ];
            dy_des_unscaled = [ swingX.vel;                     swingZ.vel;     comZ.vel;  qPelvis.vel ];
            
            y_des = [(swingX.pos+obj.stanceToePos(1)).*OutputScaling.SwingToeXZ(1);
                swingZ.pos.*OutputScaling.SwingToeXZ(2);
                comZ.pos.*OutputScaling.comXYZ(3);
                qPelvis.pos.*OutputScaling.pelvisPitch];
            dy_des = [swingX.vel.*OutputScaling.SwingToeXZ(1);
                swingZ.vel.*OutputScaling.SwingToeXZ(2);
                comZ.vel.*OutputScaling.comXYZ(3);
                qPelvis.vel.*OutputScaling.pelvisPitch];
            
            y_act_unscaled = [SwingFootPos.PosOri ;  yCOMpos.PosOri(3); yPelvisEuler.PosOri ];
            dy_act_unscaled =[SwingFootPos.VelOri ;  yCOMpos.VelOri(3); yPelvisEuler.VelOri];
            
            y_act =  [SwingFootPos.Pos;  yCOMpos.Pos(3);  yPelvisEuler.Pos];
            dy_act = [SwingFootPos.Vel;  yCOMpos.Vel(3);  yPelvisEuler.Vel];
            
            J_y =  [SwingFootPos.J;   yCOMpos.J(3,:);    yPelvisEuler.J];
            dJ_y = [SwingFootPos.dJ;  yCOMpos.dJ(3,:);   yPelvisEuler.dJ];
            
            %%%%
            y =  y_act - y_des;
            dy = dy_act - dy_des;
            
            y_unscaled = y_act_unscaled - y_des_unscaled;
            dy_unscaled = dy_act_unscaled - dy_des_unscaled;
            obj.RD2 = struct('y',  y,  'dy',  dy,  'J_y', J_y, 'dJ_y', dJ_y, 'y_des', y_des, 'dy_des',  dy_des,... %'ddy_des',ddy_des, ...
                'dy_unscaled', dy_unscaled, 'y_des_unscaled',y_des_unscaled,  'dy_des_unscaled',dy_des_unscaled, 'y_unscaled',y_unscaled);
            
            %%% DMC use next output
            N = round(obj.tSSP/obj.period) - tIdx + 2;
            % stepLnext = obj.swingXnext([SwingFootPos.PosOri(1); SwingFootPos.VelOri(1)], ...
            %                           [stepLengthDes + obj.stanceToePos(1); stepLengthVelDes], N);
            Ldes  = (stepLengthDes + obj.stanceToePos(1))*obj.smfVec(tIdx) + SwingFootPos.PosOri(1)*(1 - obj.smfVec(tIdx));
            % Ldes = stepLengthDes + obj.stanceToePos(1);
            % stepLnext = obj.swingXnext([SwingFootPos.PosOri(1); 0], ...
            %                                      [Ldes; stepLengthVelDes], N);
            
            %            stepLnext = obj.swingXnext([SwingFootPos.PosOri(1); SwingFootPos.VelOri(1)], ...
            %                                      [swingX.pos+obj.stanceToePos(1); stepLengthVelDes], N);
            
            obj.y_kp1 = [swingX.pos+obj.stanceToePos(1);% swingX.posPlus; %stepLnext(1); %obj.swingZtraj.pos(tIdx);% swingX.posPlus;
                obj.swingZtraj.pos(tIdx+1);
                obj.massHeight;
                obj.qPelvis0];
        end
        
        function obj = getCOM(obj, q, dq, ddq)
            J = J_COMPosition(q);
            dJ = dJ_COMPosition(q, dq);
            accCOM = J*ddq + dJ*dq;
            velCOM = J*dq;
            obj.accCOMvec(1:end-1) = obj.accCOMvec(2:end);
            obj.velCOMvec(1:end-1) = obj.velCOMvec(2:end);
            obj.accCOMvec(end) = accCOM(1);
            obj.velCOMvec(end) = velCOM(1);
        end
        
        function [stepLengthDes, stepLengthVelDes] = LIPbasedStepping(obj, q, dq, Contact)
            %%%
            com = COMPosition(q);
            comvel = COM_velocity(q, dq);
            p = com(1) - obj.stanceToePos(1);
            v = comvel(1);
            %v = mean(obj.velCOMvec);
            a = mean(obj.accCOMvec);
            xnow = [p;v];
            K = [1.2444, 0.3498];
            K = obj.Kdeadbeat;
            %%%% based on the H-lip plan step size
            switch obj.gaitType
                case 'P1'
                    S = transpose((obj.Ahat0 - eye(2))*obj.Xdes)*[1;0]/obj.Bhat0(1);
                    % S = vxdes/([0,1]*inv(obj.Anow - eye(2))*obj.Bnow);  %%% the regressed formulation %%%%%%%%%
                    %%%% in the form for feedback gain*error
                    stepLengthDes = K*(xnow - obj.Xdes) - S;
                case 'P2'
                    if  Contact == AmberConstants.LeftFootContact
                        xdes = obj.XLdes;
                    else
                        xdes = obj.XRdes;
                    end
                    TD = 0;
                    stepLengthNominal = 2*xdes(1) + TD*xdes(2);
                    stepLengthDes = K*(xnow - xdes) + stepLengthNominal;
            end
            stepLengthVelDes =  K*[v;a];
            
            obj.pVec = [obj.pVec, p];
            obj.vVec = [obj.vVec, v];
            obj.uVec = [obj.uVec, stepLengthDes]; 
        end
        
        function [swingX, stepLengthDes, stepLengthVelDes] = swingFootPlanning(obj, q, dq, tIdx, Contact)
            
            [stepLengthDes, stepLengthVelDes] = obj.LIPbasedStepping(q, dq, Contact) ;
            %%% smoothing from the current swing foot position
%             if  Contact == AmberConstants.LeftFootContact
%                 stepLengthPre = pRightToe(q) - pLeftToe(q);
%                 stepLengthPreVel = vRightToe([q;dq]);
%             elseif  Contact == AmberConstants.RightFootContact
%                 stepLengthPre = pLeftToe(q) - pRightToe(q);
%                 stepLengthPreVel = vLeftToe([q;dq]);
%             end
%             stepLengthPre = stepLengthPre(1);
%             stepLengthPreVel = stepLengthPreVel(1);
            
            cNow  = obj.smfVec(tIdx); %cNow from 0 ->1
            dcNow = obj.smfVecVel(tIdx);
            c2Now = obj.smfVec2(tIdx);
            intc2 = cumtrapz(obj.smfVec2)*obj.period;
            dc2Now = obj.smfVecVel2(tIdx);
            
            p0 = obj.SSPnsf0(1) - obj.stanceToePos(1);
            v0 = obj.SSPnsfVel0(1);
            p_down = p0 + v0*tIdx*obj.period - v0*intc2(tIdx);
            v_down = v0*(1 - c2Now);
            stepL = stepLengthDes*cNow + (1-cNow)*p_down ; % (1-c2Now)*stepLengthPreVel*tIdx*obj.period ;
            stepLvel = stepLengthVelDes*cNow + stepLengthDes*dcNow - dcNow*p_down +  (1-cNow)*v_down; % +...
            %  (1-c2Now)*stepLengthPreVel + (-dc2Now)*stepLengthPreVel*tIdx*obj.period ;
            %  stepLvel = stepLengthDes*dcNow - dcNow*stepLengthPre;
            
            cNext = obj.smfVec(tIdx);
            p_down_next = p0 + v0*(tIdx+1)*obj.period - v0*intc2(tIdx+1);
            stepLnext =  stepLengthDes*cNext + (1-cNext)*p_down_next;
            swingX = struct('pos', stepL, 'vel', stepLvel, 'posPlus', stepLnext);
        end
        
        function out = swingXnext(obj, x1, xTarget, N)
            h = obj.period;
            % X1 is the current state
            % Xtarget is the target
            
            At = [1, h; 0, 1];
            Bt = [h^2/2; h];
            
            %%% try QP %%%%
            Xvar = sdpvar(2,   N);
            Uvar = sdpvar(N-1, 1);
            constraints = [Xvar(:, 1) == x1;
                Xvar(1, N) == xTarget(1)];
            
            for i = 1:N-1
                constraints = [constraints; Xvar(:, i+1) == At*Xvar(:, i) + Bt*Uvar(i)];
            end
            % cost = (Xvar(1,:) - xTarget(1))*(Xvar(1,:) - xTarget(1))' +  Uvar'*Uvar;
            % cost =  Uvar'*Uvar;
            cost =  Xvar(2,:)*Xvar(2,:)'; %+ 0*Uvar'*Uvar;
            optimize(constraints, cost);
            out = double(Xvar(:,2));
            
            Xvar = double(Xvar);
            figure(1234321), box on, grid on;
            plot(1:N, Xvar(1,:), 'r-o'), ylim([-1,1]); xlim([0, 40]);
        end 
            
        function ykp1 = nextDesiredOutput(obj, q, dq, tIdx, Contact) 
             %  swingX = obj.swingFootPlanning(q, dq, tIdx, Contact);
              ykp1 = [ obj.swingX.posPlus; % obj.swingZtraj.pos(tIdx+1); %swingX.pos; %+obj.stanceToePos(1)) ;
                obj.swingZtraj.pos(tIdx+1);
                obj.massHeight;
                obj.qPelvis0];
        end 
        
        function SwingFootZconstruct(obj)
            %%% create swing foot trajectory
            %%% stepX, stepY: step size
            Tf = obj.tSSP;
            Tvec = 0:obj.period:Tf;
            HeightZ = obj.swingHeight;
            %%%% swingZ, velocity is a sinusoidal
            swingZvel = HeightZ*pi/Tf*sin(Tvec/Tf*2*pi);
            makeSureIntrusion = 0.005;
            swingZ = HeightZ/2*(-cos(Tvec/Tf*2*pi)+1) - makeSureIntrusion*Tvec/Tf;
            
            vimpact = obj.impactVel; %0.05;
            %%%%% add velocity adjustment
            N = length(Tvec);
            TvecQuad1 = Tvec(1:floor(N/4));
            TvecQuad2 = Tvec(ceil(N/4):floor(N/2));
            TvecHalf2 = Tvec(ceil(N/2):N);
            
            vel_1 = vimpact/(Tf/4)*TvecQuad1;
            vel_2 = vimpact - vimpact/(Tf/4)*(TvecQuad2 - Tf/4);
            vel_3 = -vimpact/(Tf/2)*(TvecHalf2-Tf/2);
            addZvel =  -[vel_1, vel_2, vel_3];
            addZvel = addZvel(1:N);
            
            swingZvel = swingZvel + addZvel;
            swingZpos =  swingZ + cumsum(addZvel)*obj.period;
            
            addIntrusionVel = ones(1,100)*swingZvel(end);
            
            swingZpos = [swingZpos,  swingZpos(end)+cumsum(addIntrusionVel)*obj.period];
            swingZvel = [swingZvel, addIntrusionVel];
            obj.swingZtraj = struct('pos', swingZpos, 'vel', swingZvel);
        end
        
    end
    
    
    methods %%%% HZD type
        % function VRreconstruct(obj)
        %     %%% however vr is based on time:::
        %     a = obj.a;
        %     T = obj.p(1)*2;
        % 
        %     s = 0:0.001:1;
        %     tvec = linspace(0, T, 2*length(s));
        %     pos = [];
        %     vel = [];
        %     for i =1:length(s)
        %         pos = [pos, bezier(a, s(i))];
        %         vel = [vel, dbezier(a, s(i))];
        %     end
        %     pos_switch_stance = [pos(3,:); pos(4,:); pos(1,:); pos(2,:)];
        %     vel_switch_stance = [vel(3,:); vel(4,:); vel(1,:); vel(2,:)];
        % 
        %     obj.tvec = tvec;
        %     obj.qmdes = [pos, pos_switch_stance];
        %     obj.dqmdes = [vel, vel_switch_stance];
        % end 
        % VIRTUAL CONSTRAINT
        
        function timeBasedConstruct(obj, logOpt)
            tvec = logOpt.t;
            
            pos = logOpt.qR(4:7, :);
            vel = logOpt.dqR(4:7, :);
            
            pos_switch_stance = [pos(3,:); pos(4,:); pos(1,:); pos(2,:)];
            vel_switch_stance = [vel(3,:); vel(4,:); vel(1,:); vel(2,:)];
            
            obj.tvec = [tvec, tvec(end) + tvec(2:end)];
            obj.qmdes = [pos, pos_switch_stance(:, 2:end)];
            obj.dqmdes = [vel, vel_switch_stance(:, 2:end)];
        end
        
        function VRangleConstruct(obj, logOpt)
            
            tvec = logOpt.t;
            qvec = logOpt.qR;
            dqvec = logOpt.dqR;
            
            N = length(tvec);
            hipPosVec = zeros(1, N);
            hipVelvec = zeros(1, N);
            for i = 1:length(tvec)
                qi = qvec(:, i);
                dqi = dqvec(:, i);
                hipPos = pHip(qi);
                hipVel = vHip([qi; dqi]);
                
                hipPosVec(i) = hipPos(1);
                hipVelVec(i) = hipVel(1);
            end
            
            obj.p0 = hipPosVec(1);
            obj.pN = hipPosVec(end);
            
            sNormalized = (hipPosVec -  obj.p0)/(obj.pN -  obj.p0);
            qm = qvec(4:7,:);
            dqm = dqvec(4:7, :);
            ndeg = 6;
            pmat = zeros(4, ndeg+1);
            for i = 1:4
                pmat(i, :) = polyfit(sNormalized, qm(i, :), ndeg);
            end
            obj.polyMat = pmat;
            
            
        end
        
        function getDesiredOutputs(obj, q, dq, t, contact)
            switch obj.timeOrPhase
                case 'time'
                    obj.getDesiredOutputsTime(q, dq, t, contact);
                case 'phase'
                    obj.getDesiredOutputsPhase(q, dq, t, contact);
            end
        end
        
        function getDesiredOutputsTime(obj, q, dq, t, contact)
            
            %             s = mod(t, (obj.p(1) - obj.p(2)));
            %             s = (s- obj.p(2))/ (obj.p(1) - obj.p(2));
            %
            % %             M = 5; %  alpha = sym('alpha%d%d', [4 M+1]);
            % %             Bezier_bs = @(s,alpha,i,k) alpha(i,k+1)*factorial(M)/(factorial(k)*factorial(M-k))*(s^k)*((1-s)^(M-k));
            % %             Bezier_qN_qq = @(s,alpha,i) 0;
            % %             for j = 0:M
            % %                 Bezier_qN_qq = @(s,alpha,i) Bezier_bs(s,alpha,i,j) + Bezier_qN_qq(s,alpha,i);
            % %             end
            % %             Bezier_qN = @(s,alpha) [Bezier_qN_qq(s,alpha,1);Bezier_qN_qq(s,alpha,2);...
            % %                 Bezier_qN_qq(s,alpha,3);Bezier_qN_qq(s,alpha,4)];
            %
            % %             out1 =   Bezier_qN_qq(s,alpha,1);
            % %             out2 =   Bezier_qN_qq(s,alpha,2);
            % %             out3 =   Bezier_qN_qq(s,alpha,3);
            % %             out4 =   Bezier_qN_qq(s,alpha,4);
            % %
            %             qmdes = bezier(obj.a, s);
            %             dqmdes = dbezier(obj.a, s);
            %             if ~strcmp(contact, AmberConstants.RightFootContact)
            %                 qmdes = [qmdes(3); qmdes(4); qmdes(1); qmdes(2)];
            %                 dqmdes = [dqmdes(3); dqmdes(4); dqmdes(1); dqmdes(2)];
            %             end
            
            t = mod(t, obj.tvec(end));
            qmdes = zeros(4,1); dqmdes = zeros(4,1);
            for i = 1:4
                qmdes(i) = interp1(obj.tvec, obj.qmdes(i,:), t);
                dqmdes(i) = interp1(obj.tvec, obj.dqmdes(i,:), t);
            end
            
            y_des_unscaled =  qmdes;
            dy_des_unscaled = dqmdes;
            
            y_des = y_des_unscaled;
            dy_des = dy_des_unscaled;
            
            y_act_unscaled = q(obj.motorIdx);
            dy_act_unscaled = dq(obj.motorIdx);
            
            y_act = y_act_unscaled;
            dy_act = dy_act_unscaled;
            
            J_y =  obj.Jmotors;
            dJ_y = obj.dJmotors;
            
            %%%%
            y =  y_act - y_des;
            dy = dy_act - dy_des;
            
            y_unscaled = y_act_unscaled - y_des_unscaled;
            dy_unscaled = dy_act_unscaled - dy_des_unscaled;
            obj.RD2 = struct('y',y,'dy',dy,'J_y', J_y, 'dJ_y', dJ_y,'y_des',y_des,'dy_des',dy_des,... %'ddy_des',ddy_des, ...
                'dy_unscaled', dy_unscaled, 'y_des_unscaled',y_des_unscaled,  'dy_des_unscaled',dy_des_unscaled, 'y_unscaled',y_unscaled);
            
        end
        
        function getDesiredOutputsPhase(obj, q, dq, t, contact)
            
            switch contact
                case AmberConstants.RightFootContact
                    pstance = pRightToe(q);
                case AmberConstants.LeftFootContact
                    pstance = pLeftToe(q);
            end
            hipPos = pHip(q) - pstance;
            hipVel = vHip([q; dq]);
            qmdes = zeros(4,1); dqmdes = zeros(4,1);
            s = (hipPos(1) - obj.p0)/(obj.pN - obj.p0);
            ds = 1/(obj.pN - obj.p0)*hipVel(1);
            
            for i = 1:4
                qmdes(i) = polyval( obj.polyMat(i, :), s);
                dpcoef = polyder(obj.polyMat(i, :));
                dqmdes(i) = polyval( dpcoef, s).*ds;
            end
            
            if  contact == AmberConstants.LeftFootContact
                qmdes = [qmdes(3); qmdes(4); qmdes(1); qmdes(2)];
                dqmdes = [dqmdes(3); dqmdes(4); dqmdes(1); dqmdes(2)];
            end
            
            y_des_unscaled =  qmdes;
            dy_des_unscaled = dqmdes;
            
            y_des = y_des_unscaled;
            dy_des = dy_des_unscaled;
            
            y_act_unscaled = q(obj.motorIdx);
            dy_act_unscaled = dq(obj.motorIdx);
            
            y_act = y_act_unscaled;
            dy_act = dy_act_unscaled;
            
            J_y =  obj.Jmotors;
            dJ_y = obj.dJmotors;
            
            %%%%
            y =  y_act - y_des;
            dy = dy_act - dy_des;
            
            y_unscaled = y_act_unscaled - y_des_unscaled;
            dy_unscaled = dy_act_unscaled - dy_des_unscaled;
            obj.RD2 = struct('y',y,'dy',dy,'J_y', J_y, 'dJ_y', dJ_y,'y_des',y_des,'dy_des',dy_des,... %'ddy_des',ddy_des, ...
                'dy_unscaled', dy_unscaled, 'y_des_unscaled',y_des_unscaled,  'dy_des_unscaled',dy_des_unscaled, 'y_unscaled',y_unscaled);
            
        end
    end
    
    methods % output functions
        
        function outputSwingFootPos = OutputSwingFootPos(obj, q, dq, OutputScaling, Contact)
            % of course it should be rigid foot pos
            x = [q;dq];
            outputSwingFootPos = struct('PosOri', [], 'VelOri', [], 'Pos', [], 'Vel', [], 'J', [], 'dJ', []);
            if  Contact == AmberConstants.LeftFootContact
                % scale SwingFoot Pos
                SwingFootPosUnscaled = pRightToe(q);
                SwingFootVelUnscaled = vRightToe(x);
                SwingFootPos = SwingFootPosUnscaled.*OutputScaling.SwingToeXZ;
                SwingFootVel = SwingFootVelUnscaled.*OutputScaling.SwingToeXZ;
                J_SwingFootPos  = J_rightToe(q).*OutputScaling.SwingToeXZ;
                dJ_SwingFootPos = dJ_rightToe(x).*OutputScaling.SwingToeXZ;
            else
                % scale SwingFoot Pos
                SwingFootPosUnscaled = pLeftToe(q);
                SwingFootVelUnscaled = vLeftToe(x);
                SwingFootPos = SwingFootPosUnscaled.*OutputScaling.SwingToeXZ;
                SwingFootVel = SwingFootVelUnscaled.*OutputScaling.SwingToeXZ;
                J_SwingFootPos  = J_leftToe(q).*OutputScaling.SwingToeXZ;
                dJ_SwingFootPos = dJ_leftToe(x).*OutputScaling.SwingToeXZ;
            end
            outputSwingFootPos.PosOri = SwingFootPosUnscaled;
            outputSwingFootPos.VelOri = SwingFootVelUnscaled;
            outputSwingFootPos.Pos = SwingFootPos;
            outputSwingFootPos.Vel = SwingFootVel;
            outputSwingFootPos.J = J_SwingFootPos;
            outputSwingFootPos.dJ = dJ_SwingFootPos;
        end
        
        function outputPelvisEuler = OutputPelvisEuler(obj, q, dq, OutputScaling)
            outputPelvisEuler = struct('PosOri', [], 'VelOri', [], 'Pos', [], 'Vel', [], 'J', [], 'dJ', []);
            % scale pelvis orientation
            Pelvis_ori_unscaled = pelvis_ori(q);  %3x1
            Pelvis_rotV_unscaled = pelvis_rotvel(q,dq);
            
            % scale pelvis ori
            Pelvis_Ori = Pelvis_ori_unscaled.*OutputScaling.pelvisPitch;
            Pelvis_rotV = Pelvis_rotV_unscaled.*OutputScaling.pelvisPitch;
            J_Pelvis_Ori = J_pelvis_ori(q).*OutputScaling.pelvisPitch;
            dJ_Pelvis_Ori = dJ_pelvis_ori(q,dq).*OutputScaling.pelvisPitch;
            
            outputPelvisEuler.PosOri = Pelvis_ori_unscaled;
            outputPelvisEuler.VelOri = Pelvis_rotV_unscaled;
            outputPelvisEuler.Pos = Pelvis_Ori;
            outputPelvisEuler.Vel = Pelvis_rotV;
            outputPelvisEuler.J = J_Pelvis_Ori;
            outputPelvisEuler.dJ = dJ_Pelvis_Ori;
        end
        
        function outputCOM = OutputCOMpos(obj,q ,dq ,OutputScaling)
            %%% relative to the stance toe
            x = [q;dq];
            outputCOM = struct('PosOri', [], 'VelOri', [], 'Pos', [], 'Vel', [], 'J', [], 'dJ', []);
            %             if strcmp(Contact, AmberConstants.LeftFootContact)
            %                 % scale SwingFoot Pos
            %                 StanceFootPos = pLeftToe(q)'; % rigidFootPos(q)'*[0;1];
            %             elseif strcmp(Contact, AmberConstants.RightFootContact)
            %                 % scale SwingFoot Pos
            %                 StanceFootPos = pRightToe(q)';
            %             end
            posUnscaled = COMPosition(q); % - StanceFootPos;
            velUnscaled = COM_velocity(q,dq);
            comPos = posUnscaled.*OutputScaling.comXYZ';
            comVel = velUnscaled.*OutputScaling.comXYZ';
            J_com  = J_COMPosition(q).*OutputScaling.comXYZ;
            dJ_com = dJ_COMPosition(q, dq).*OutputScaling.comXYZ;
            outputCOM.PosOri = posUnscaled;
            outputCOM.VelOri = velUnscaled;
            outputCOM.Pos = comPos;
            outputCOM.Vel = comVel;
            outputCOM.J = J_com;
            outputCOM.dJ = dJ_com;
        end
        
        function [y_sol,yd_sol,ydd_sol] = LIP_sol(obj, lambda,y0,yd0,t)
            
            %%%%%%%%%%%% y_ddt = g/z0 * y   % LIP dynamics
            %%% y(t) = c1 e^(v*t) + c2 e^(-v*t)
            %%% y_dot(t) = v( c1 e^(v*t) - c2 e^(-v*t))
            
            c1 = @(y0,yd0) 1/2*(y0 + yd0/lambda);
            c2 = @(y0,yd0) 1/2*(y0 - yd0/lambda);
            
            f_y_sol  = @(t, y0, yd0) c1(y0, yd0)*exp(lambda*t) + c2(y0, yd0)*exp(-lambda*t);
            f_yd_sol = @(t, y0, yd0) lambda*(c1(y0, yd0)*exp(lambda*t) - c2(y0, yd0)*exp(-lambda*t));
            f_ydd_sol = @(t, y0, yd0) lambda^2 * y0;
            if length(t) == 1
                y_sol = f_y_sol(t,y0, yd0);
                yd_sol = f_yd_sol(t,y0, yd0);
                ydd_sol = f_ydd_sol(t,y0, yd0);
            else
                y_sol = zeros(size(t));
                yd_sol = zeros(size(t));
                ydd_sol = zeros(size(t));
                
                for i = 1:length(t)
                    y_sol(i) =  f_y_sol(t(i),y0, yd0);
                    yd_sol(i) =  f_yd_sol(t(i),y0, yd0);
                    ydd_sol(i) =  f_ydd_sol(t(i),y0, yd0);
                end
            end
        end
        
    end
    
    methods %%% past
        % function swingX = swingFootPlanningP1(obj, q, dq, tIdx, Contact)
        %     %%%
        %     com = COMPosition(q);
        %     comvel = COM_velocity(q, dq);
        %     p = com(1) - obj.stanceToePos(1);
        %     v = comvel(1);
        %     a = obj.g/obj.massHeight*p;
        %     %v = mean(obj.velCOMvec);
        %     a = mean(obj.accCOMvec);
        %     xnow = [p;v];
        %     %%%% based on the H-lip plan step size
        %     S = transpose((obj.Ahat0 - eye(2))*obj.Xdes)*[1;0]/obj.Bhat0(1);
        %     % S = vxdes/([0,1]*inv(obj.Anow - eye(2))*obj.Bnow);  %%% the regressed formulation %%%%%%%%%
        %     %%%% in the form for feedback gain*error
        %     K = [1.2444, 0.3498];
        %     K = obj.Kdeadbeat;
        % 
        %     L = K*(xnow - obj.Xdes);
        %     stepLengthDes = L - S;
        %     stepLengthVelDes =  K*[v;a];
        %     %%% smoothing from the current swing foot position
        %     if  Contact == AmberConstants.LeftFootContact
        %         stepLengthPre = pRightToe(q) - pLeftToe(q);
        %         stepLengthPreVel = vRightToe([q;dq]);
        %     elseif  Contact == AmberConstants.RightFootContact
        %         stepLengthPre = pLeftToe(q) - pRightToe(q);
        %         stepLengthPreVel = vLeftToe([q;dq]);
        %     end
        %     stepLengthPre = stepLengthPre(1);
        %     stepLengthPreVel = stepLengthPreVel(1);
        %     %             stepLengthPre = obj.SSPnsf0(1) - obj.stanceToePos(1);
        %     %             stepLengthPreVel = obj.SSPnsfVel0(1);
        %     cNow = obj.smfVec(tIdx); %cNow from 0 ->1
        %     dcNow = obj.smfVecVel(tIdx);
        %     c2Now = obj.smfVec2(tIdx);
        %     intc2 = cumtrapz(obj.smfVec2)*obj.period;
        %     dc2Now = obj.smfVecVel2(tIdx);
        %     stepL = stepLengthDes*cNow + (1-cNow)*stepLengthPre ; % (1-c2Now)*stepLengthPreVel*tIdx*obj.period ;
        %     stepLvel = stepLengthVelDes*cNow + stepLengthDes*dcNow - dcNow*stepLengthPre +  (1-cNow)*stepLengthPreVel; % +...
        %     %  (1-c2Now)*stepLengthPreVel + (-dc2Now)*stepLengthPreVel*tIdx*obj.period ;
        %     %  stepLvel = stepLengthDes*dcNow - dcNow*stepLengthPre;
        % 
        %     p0 = obj.SSPnsf0(1) - obj.stanceToePos(1);
        %     v0 = obj.SSPnsfVel0(1);
        %     p_down = p0 + v0*tIdx*obj.period - v0*intc2(tIdx);
        %     v_down = v0*(1 - c2Now);
        %     stepL = stepLengthDes*cNow + (1-cNow)*p_down ; % (1-c2Now)*stepLengthPreVel*tIdx*obj.period ;
        %     stepLvel = stepLengthVelDes*cNow + stepLengthDes*dcNow - dcNow*p_down +  (1-cNow)*v_down; % +...
        %     %  (1-c2Now)*stepLengthPreVel + (-dc2Now)*stepLengthPreVel*tIdx*obj.period ;
        %     %  stepLvel = stepLengthDes*dcNow - dcNow*stepLengthPre;
        % 
        %     swingX = struct('pos', stepL, 'vel', stepLvel);
        % 
        %     obj.pVec = [obj.pVec, p];
        %     obj.vVec = [obj.vVec, v];
        % end
        
    end
end

