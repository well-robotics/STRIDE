classdef analysisWalker < handle 
    %%%% generation plots for analysis
    properties 
        log 
        flow
        outputs
        Ncycles
        eventSim 
        freq
        
    end 
    
    methods 
        function obj = analysisWalker(log, outputs, Ncycles, eventSim, freq) 
            obj.log = log; 
            obj.flow = log.flow; 
            obj.outputs = outputs;
            obj.Ncycles = Ncycles;
            % obj.eventSim = eventSim;
            obj.initialize();
            obj.freq = freq;
            
        end 
        
        function initialize(obj)
            global colorString MotorList  floatingBase motorIdx baseIdx
            colorString = {'r','g','b','k','m','c','y','w'};
            MotorList = {' LeftHipPitch ',' LeftKneePitch ',' RightHipPitch ',' RightKneePitch '};
            floatingBase = {' PelvisX ',' PelvisZ ',' Pelvis Pitch '};
            motorIdx = 4:7;
            baseIdx = 1:3;
           
            set(0, 'DefaultLineLineWidth', 1);
        end
        
        function animate(obj) 
            log = obj.log;            
            t = log.flow.t;
            q = log.flow.q;
            
            anim = Animator.AmberAnimator(t, q');
            anim.pov = Animator.AnimatorPointOfView.West;
            anim.Animate(true);
            anim.isLooping = false;
            anim.updateWorldPosition = true;
            
            conGUI = Animator.AnimatorControls();
            conGUI.anim = anim;
           
          %  anim.axisOption = 'fixed';
           % anim.timeTextOption = 'fixed';
           % anim.timeTextLocation = [-3, 3, 1.2];
            % anim.axisSize = [-5, 1.5, -4, 4,-0.5, 1.5];
           
        end 
    end 
    
    methods %%% MPC shit
        % function visualize(obj)
        %     log = obj.log;
        %     t = log.flow.t;
        %     q = log.flow.q;
        % 
        %     N = length(t); 
        %     pLknee = zeros(N,3); 
        %     pRknee = zeros(N,3);
        %     pLtoe = zeros(N,3); 
        %     pRtoe = zeros(N,3);
        %     phip = zeros(N,3); 
        %     for i = 1:N
        %         pLknee(i,:) = pLeftKnee(q(i,:)); 
        %         pRknee(i,:) = pRightKnee(q(i,:)); 
        %         pLtoe(i,:) = pLeftToe(q(i,:));
        %         pRtoe(i,:) = pRightToe(q(i,:));
        %         phip(i,:) = pHip(q(i,:)); 
        %     end
        %     data = struct; 
        %     data.t = {t};
        %     data.mx = phip(:,1); 
        %     data.mz = phip(:,3); 
        %     data.qT = q(:,3); 
        %     data.qRt = q(:,4); 
        %     data.qRs = q(:,5); 
        %     data.qLt = q(:,6); 
        %     data.qLs = q(:,7); 
        %     data.leftKneeX = pLknee(:,1); 
        %     data.rightKneeX = pRknee(:,1); 
        %     data.leftKneeZ = pLknee(:,3); 
        %     data.rightKneeZ = pRknee(:,3); 
        % 
        %     data.fLx = pLtoe(:,1); 
        %     data.fLz = pLtoe(:,3); 
        %     data.fRx = pRtoe(:,1); 
        %     data.fRz = pRtoe(:,3); 
        % 
        %     stones = struct; 
        %     stones.width = 0.1;
        %     stones.height = [0.1, 0.2, 0.1];
        %     stones.distance = [0.2, 0.4, 0.1];
        %     scene = FiveLinkWalkerScene(data, stones);
        %     Player(scene);
        % end
        % 
        % function checkMPCsol(obj, mpc, target)
        %     sol = mpc.solLog; 
        %     nx = (mpc.N+1)*mpc.nx; 
        %     nu = mpc.nu*mpc.N; 
        %     nF = mpc.nF*mpc.N; 
        %     xsol = sol(1:nx, :); 
        %     usol = sol(nx+1:nx+nu, :); 
        %     Fsol = sol(nx+nu+1:end-1, :); 
        %     tvec = sol(end,:); 
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     sim = obj.log.flow; 
        %     des = target.flow;
        %     figure, 
        %     for i = 1:7
        %         subplot(2, 7, i)
        %         plot(tvec, xsol(i,:), 'r', des.t, des.q(:,i), 'b', sim.t, sim.q(:,i), 'g');
        %         xlim([0, sim.t(end)]); 
        %     end
        %     legend('mpc','nominal', 'sim');
        %     for i = 1:7
        %         subplot(2, 7, 7+i)
        %         plot(tvec, xsol(7+i,:), 'r', des.t, des.dq(:,i), 'b', sim.t, sim.dq(:,i), 'g');    
        %         xlim([0, sim.t(end)]);
        %     end
        %     title('states');
        % 
        %     %%%%%%%%%%%%%%% torque comp %%%%%%%%%%%%%%%%%%
        %     figure,
        %     plot(tvec, usol(1:4,:), '-', des.t, des.u, '--');
        %     title('torque');
        %     legend('mpc sol', 'nominal');
        %     xlim([0, sim.t(end)]);
        % 
        %     %%%%%%%%%%%%%%% GRF comp %%%%%%%%%%%%%%%%%%%%%
        %     figure,
        %     subplot(2, 1, 1);
        %     plot(tvec, Fsol(1,:), '-', sim.t, sim.F(1,:), '--');
        %     title('forward force'); legend('mpc sol', 'simulation');
        %     subplot(2, 1, 2);
        %     plot(tvec, Fsol(2,:), '-', sim.t, sim.F(2,:), '--');
        %     title('normal force'); legend('mpc sol', 'simulation');
        %     xlim([0, sim.t(end)]);
        % 
        %     figure,
        %     subplot(2, 1, 1);
        %     plot(tvec, Fsol(1,:), 'r-', tvec, Fsol(3,:), 'b-',...
        %         des.t, target.sim.leftFootForce(1,:), 'r--', des.t, target.sim.rightFootForce(1,:), 'r--');
        %     title('forward force'); legend('mpc sol', 'desired');
        %     xlim([0, sim.t(end)]);
        % 
        %     subplot(2, 1, 2);
        %     plot(tvec, Fsol(2,:), 'r-', tvec, Fsol(4,:), 'b-', ...
        %         des.t, target.sim.leftFootForce(2,:), 'r--', des.t, target.sim.rightFootForce(2,:), 'b--');
        %     title('normal force'); legend('mpc sol', 'desired');
        %     xlim([0, sim.t(end)]);
        % end
        % 
        % function checkMPCsolHorizon(obj, mpc, target, t)
        %     idx = round(t*1000/mpc.Nspace)+1; 
        %     sol = mpc.solLog; 
        %     nx = (mpc.N+1)*mpc.nx; 
        %     nu = mpc.nu*mpc.N; 
        %     nF = mpc.nF*mpc.N; 
        %     xsol = sol(1:nx, idx); 
        %     usol = sol(nx+1:nx+nu, idx); 
        %     Fsol = sol(nx+nu+1:end-1, idx); 
        %     tIdx = sol(end,idx) + 0.001; 
        %     tvec = tIdx:0.001*mpc.Nspace:tIdx+(mpc.N-1)*0.001*mpc.Nspace;
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     sim = obj.log.flow; 
        %     des = target.flow;
        %     xsol = reshape(xsol, mpc.nx, []); 
        %     usol = reshape(usol, mpc.nu, []); 
        %     Fsol = reshape(Fsol, mpc.nF, []); 
        %     figure, 
        %     for i = 1:7
        %         subplot(2, 7, i)
        %         plot(tvec, xsol(i,1:end-1), 'r', des.t, des.q(:,i), 'b', sim.t, sim.q(:,i), 'g');
        %         xlim([tvec(1), tvec(end)]); 
        %     end
        %     legend('mpc','nominal', 'sim');
        %     for i = 1:7
        %         subplot(2, 7, 7+i)
        %         plot(tvec, xsol(7+i,1:end-1), 'r', des.t, des.dq(:,i), 'b', sim.t, sim.dq(:,i), 'g');    
        %         xlim([tvec(1), tvec(end)]); 
        %     end
        %     title('states');
        % 
        %     %%%%%%%%%%%%%%% torque comp %%%%%%%%%%%%%%%%%%
        %     figure,
        %     for i = 1:4
        %         subplot(4,1,i)
        %         plot(tvec, usol(i,:), '-', des.t, des.u(:,i), '--');
        %         title('torque');
        %         legend('mpc sol', 'nominal');
        %         xlim([tvec(1), tvec(end)]);
        %     end
        % 
        %     %%%%%%%%%%%%%%% GRF comp %%%%%%%%%%%%%%%%%%%%%
        %     figure,
        %     subplot(2, 1, 1);
        %     plot(tvec, Fsol(1:2:3,:), '-', sim.t, sim.F(1,:), '--');
        %     title('forward force'); legend('mpc sol', 'simulation');   
        %     xlim([tvec(1), tvec(end)]); 
        % 
        %     subplot(2, 1, 2);
        %     plot(tvec, Fsol(2:2:4,:), '-', sim.t, sim.F(2,:), '--');
        %     title('normal force'); legend('mpc sol', 'simulation');
        %     xlim([tvec(1), tvec(end)]);
        % 
        %     figure,
        %     subplot(2, 1, 1);
        %     plot(tvec, Fsol(1,:), 'r-', tvec, Fsol(3,:), 'b-',...
        %         des.t, mpc.nominalForces(1,:), 'r--', des.t, mpc.nominalForces(3,:), 'b--');
        %     title('forward force'); legend('mpc sol', 'ref');
        %     xlim([tvec(1), tvec(end)]);
        % 
        %     subplot(2, 1, 2);
        %     plot(tvec, Fsol(2,:), 'r-', tvec, Fsol(4,:), 'b-', ...
        %         des.t, mpc.nominalForces(2,:), 'r--', des.t, mpc.nominalForces(4,:), 'b--');
        %     title('normal force'); legend('mpc sol', 'ref');
        %     xlim([tvec(1), tvec(end)]);
        % end
        % 
        function checkContactForce(obj, des)
           tvec = des.flow.t; 
           leftForce = des.sim.leftFootForce; 
           rightForce = des.sim.rightFootForce; 
           
           Navg  = 5; 
           leftZ = movmean(leftForce(2,:), Navg); 
           leftX = movmean(leftForce(1,:), Navg); 
           
           rightZ = movmean(rightForce(2,:), Navg); 
           rightX = movmean(rightForce(1,:), Navg); 
           
           figure,
           subplot(2,1,1),plot(tvec, leftForce(2,:), 'r', tvec, rightForce(2,:), 'b', ...
                 tvec, leftZ, 'r--', tvec, rightZ, 'b--'); title('Force Z'); 
           subplot(2,1,2),plot(tvec, leftForce(1,:), 'r', tvec, rightForce(1,:), 'b', ...
                 tvec, leftX, 'r--', tvec, rightX, 'b--'); title('Force X'); 
        end
        
        % function checkMPC1sol(obj, mpc, target)
        %     sol = mpc.solLog; 
        %     nx = (mpc.N+1)*mpc.nx; 
        %     nu = mpc.nu*mpc.N; 
        %     nF = mpc.nF*mpc.N; 
        %     xsol = sol(1:nx, :); 
        %     usol = sol(nx+1:nx+nu, :); 
        %     Fsol = sol(nx+nu+1:end-1, :); 
        %     tvec = sol(end,:); 
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     sim = obj.log.flow; 
        %     des = target.flow;
        %     figure, 
        %     for i = 1:7
        %         subplot(2, 7, i)
        %         plot(tvec, xsol(i,:), 'r', des.t, des.q(:,i), 'b', sim.t, sim.q(:,i), 'g');
        %         xlim([0, sim.t(end)]); 
        %     end
        %     legend('mpc','nominal', 'sim');
        %     for i = 1:7
        %         subplot(2, 7, 7+i)
        %         plot(tvec, xsol(7+i,:), 'r', des.t, des.dq(:,i), 'b', sim.t, sim.dq(:,i), 'g');    
        %         xlim([0, sim.t(end)]);
        %     end
        %     title('states');
        % 
        %     %%%%%%%%%%%%%%% torque comp %%%%%%%%%%%%%%%%%%
        %     figure,
        %     plot(tvec, usol(1:4,:), '-', des.t, des.u, '--');
        %     title('torque');
        %     legend('mpc sol', 'nominal');
        %     xlim([0, sim.t(end)]);
        % 
        %     %%%%%%%%%%%%%%% GRF comp %%%%%%%%%%%%%%%%%%%%%
        %     figure,
        %     subplot(2, 1, 1);
        %     plot(tvec, Fsol(1,:), '-', sim.t, sim.F(1,:), '--');
        %     title('forward force'); legend('mpc sol', 'simulation');
        %     subplot(2, 1, 2);
        %     plot(tvec, Fsol(2,:), '-', sim.t, sim.F(2,:), '--');
        %     title('normal force'); legend('mpc sol', 'simulation');
        %     xlim([0, sim.t(end)]);
        % 
        %     figure,
        %     subplot(2, 1, 1);
        %     plot(tvec, Fsol(1,:), '-', des.t, target.sim.leftFootForce(1,:), 'r--', des.t, target.sim.rightFootForce(1,:), 'b--');
        %     title('forward force'); legend('mpc sol', 'desired');
        %     xlim([0, sim.t(end)]);
        % 
        %     subplot(2, 1, 2);
        %     plot(tvec, Fsol(2,:), '-', des.t, target.sim.leftFootForce(2,:), 'r--', des.t, target.sim.rightFootForce(2,:), 'b--');
        %     title('normal force'); legend('mpc sol', 'desired');
        %     xlim([0, sim.t(end)]);
        % end
        % 
        function plotKinematics(obj, target) 
            desQ = target.flow.q; 
            desDQ = target.flow.dq; 
            
            simQ = obj.log.flow.q; 
            simDQ = obj.log.flow.dq; 
            
            N = length(simQ); 
            pFootLeftDes = zeros(2,N); 
            pFootRightDes = zeros(2,N);
            pFootLeftSim = zeros(2,N); 
            pFootRightSim = zeros(2,N); 
            
            vFootLeftDes = zeros(2,N); 
            vFootRightDes = zeros(2,N);
            vFootLeftSim = zeros(2,N); 
            vFootRightSim = zeros(2,N); 
            
            for i = 1: N 
                q = simQ(i,:);  
                dq = simDQ(i,:); 
                plft = pLeftToe(q); 
                prft = pRightToe(q); 
                vlft = vLeftToe([q';dq']); 
                vrft = vRightToe([q';dq']); 
                
                pFootLeftSim(:,i) = plft(1:2:3)';
                pFootRightSim(:,i) = prft(1:2:3)';
                vFootLeftSim(:,i) = vlft(1:2:3);
                vFootRightSim(:,i) = vrft(1:2:3);
                
                q = desQ(i,:);  
                dq = desDQ(i,:);  
               
                plft = pLeftToe(q); 
                prft = pRightToe(q);  
                vlft = vLeftToe([q';dq']); 
                vrft = vRightToe([q';dq']); 
                
                pFootLeftDes(:,i) = plft(1:2:3)';
                pFootRightDes(:,i) = prft(1:2:3)';
                vFootLeftDes(:,i) = vlft(1:2:3);
                vFootRightDes(:,i) = vrft(1:2:3);
                
            end 
            
            figure, 
            tvec = (1:N)/1000;
            title('foot Z position')
            plot(tvec, pFootLeftSim(2,:), 'r', tvec, pFootLeftDes(2,:), 'r--', ... 
                    tvec, pFootRightSim(2,:), 'b', tvec, pFootRightDes(2,:), 'b--');
            legend('sim pos','ref pos')
        end 
        
    end 
    
    methods 
        
        function plotTracking(obj)
            log = obj.log;
            outputsNames = {'swingX', 'swingZ', 'comZ', 'pelvisQ'};
            %%%% leg length pelvis RPY
            figure,

            outputsAct = log.outputs.y2 + log.outputs.y2des;
            outputsVelAct = log.outputs.dy2 + log.outputs.dy2des;
            t = log.flow.t; 
            for oo = 1:4
                subplot(2, 4, oo)
                plot( ...% t, log.outputs.y2des(:,oo),'r-o',...
                    t, outputsAct(:,oo),'b-o', ...
                    t(2:end), log.outputs.yNext(1:end-1, oo), 'k-o')
                title([outputsNames{oo},' position '])
                subplot(2, 4, oo+4)

                plot(t, log.outputs.dy2des(:,oo),'r',...
                    t, outputsVelAct(:,oo),'b')
                title([outputsNames{oo}, ' velocity '])
            end
            legend('desired','actual')
        end
        
        function plotJoints(obj) 
            flow = obj.flow;
            log = obj.log; 
            global colorString MotorList PassiveJoints floatingBase motorIdx baseIdx jointIdx springIdx passiveIdx
            global N_DSP_outputs outputsLimDSP outputsLimDSPvel outputsLimSSP outputsLimSSPvel
            global motorLow motorUp
            figure,hold on, box on,grid on,
            title('Flow - MotorPos');
            for i = 1:12
                h = subplot(4,6,i); hold(h);plot(flow.t,flow.q(:,motorIdx(i)))
                plot(flow.t,ones(size(flow.t))*motorUp(i),'r',flow.t,ones(size(flow.t))*motorLow(i),'b')
                title(['flow ', MotorList{i}]);
                xlabel('Time (s)'),ylabel('rad')
            end
            title('Flow - MotorVel');
            for i = 1:12
                subplot(4,6,i+12),plot(flow.t,flow.dq(:,motorIdx(i)))
                title(['flow ',MotorList{i}]);
                xlabel('Time (s)'),ylabel('rad/s')
            end
            %
            figure,hold on, box on,grid on,
            title('QP - MotorPos');
            for i = 1:12
                subplot(2,6,i),plot(log.flow.t,log.flow.q(:,motorIdx(i)))
                title(MotorList{i});
                xlabel('Time (s)');
            end
            figure,hold on,box on, grid on,
            title('QP - MotorVel');
            for i = 1:12
                subplot(6,2,i),plot(log.flow.t,log.flow.q(:,passiveIdx(i)))
                title(PassiveJoints{i});
                xlabel('Time (s)');
            end
            % compare motor positions
            figure, hold on
            for i = 1:12
                h = subplot(2,6,i); hold(h);
                plot(log.flow.t,log.flow.q(:,motorIdx(i)),'--', flow.t,flow.q(:,motorIdx(i)),'-')
                if i == 3
                    plot(flow.t,ones(size(flow.t))*motorUp(i),'r',flow.t,ones(size(flow.t))*motorLow(i),'b')
                end
                title(MotorList{i});
                xlabel('Time (s)'),ylabel('rad')
            end
        end 

        function plotDomainDuration(obj)
            log = obj.log;
            domains = obj.domains;
            domainsVec = obj.domainsVec;
            domainDurationVec = [];
            for i = 1:length(domainsVec)
                DomainIdx = [domains{domainsVec(i)},int2str(i)];
                if isfield(log,DomainIdx)
                    if isfield(log.(DomainIdx),'T')
                        this_duration = log.(DomainIdx).T;
                    else
                        this_duration = 0;
                    end
                    domainDurationVec = [domainDurationVec, this_duration];
                end
            end
            figure, plot(domainDurationVec), title('domain duration');
        end

        function plotTorque(obj)
            log = obj.log; 
            global MotorList
            figure,
            for i = 1:size(log.flow.u,2)
                subplot(2,2,i),plot(log.flow.t,log.flow.u(:,i), 'r'), xlabel('Time (s)'), ylabel('Torque (Nm)');
                title(MotorList{i})
            end
        end
        
        function plotTorqueCompare(obj)
            l = obj.log; 
            global MotorList
            figure,
            for i = 1:size(l.flow.u,2)
                h = subplot(2,2,i);
                hold(h, 'on') 
                plot(l.flow.t,l.flow.u(:,i), 'b');
                stairs( l.flow.t,l.flow.uDM(:,i)), xlabel('Time (s)'), ylabel('Torque (Nm)');
                title(MotorList{i})
            end
        end
        
        function plotQP(obj)
            % delta CLF
            log = obj.log; 
            figure, hold on,box on
            subplot(4, 1, 1),plot(log.flow.t,log.QP.delta,'r');%ylim([0 100]);
            title('QP-delta')
            subplot(4, 1, 2),plot(log.flow.t,log.QP.exitflag,'o-')
            title('exitflag')
            subplot(4, 1, 3),plot(log.flow.t,log.QP.numiter,'g')
            title('Num of iter')
            subplot(4, 1, 4),plot(log.flow.t,log.QP.fval,'k')
            title('function value')
            
            figure, 
            subplot(2,1,1) 
            plot(log.flow.t, log.QP.GRF(1,:)), title('forward Fx')
            subplot(2,1,2) 
            plot(log.flow.t, log.QP.GRF(2,:)), title('forward Fz')
            
        end
        
         function plotDMC(obj)
            % delta CLF
            DMC = obj.log.DMC; 
            log = obj.log;
            figure, hold on, box on
            subplot(2, 2, 1), plot(log.flow.t,log.DMC.u, 'r'); 
            title('DMC-u')
            subplot(2, 2, 2),plot(log.flow.t,log.DMC.exitflag, 'o-')
            title('DMC-exitflag')
            subplot(2, 2, 3),plot(log.flow.t,log.DMC.F, 'k')
            title('DMC-contact Force')
            subplot(2, 2, 4),plot(log.flow.t,log.DMC.cost, 'k')
            title('DMC-function value')
            
          
        end
        
        % function plotGRF(obj) 
        %     eventSim = obj.eventSim;
        %     period = 0.001;
        %     %% compare foot force (midfoot force)
        %     figure,hold on,
        %     N = size(eventSim.leftFootForce,2);
        %     plot((1:N)*period, eventSim.leftFootForce(2,1:N), 'r', (1:N)*period, eventSim.rightFootForce(2,1:N), 'b');
        %     xlabel('Time (s)'),ylabel('Ground Force z(N)');
        % end
        
        function plotCOM(obj) 
            log = obj.log;
            outputs = obj.outputs;
            %             figure,
            %             subplot(2,1,1)
            %             plot(log.flow.t, log.flow.comX(1,:))
            %             title('com pos')
            %             subplot(2,1,2)
            %             plot(log.flow.t, log.flow.comX(2,:))
            %             title('com vel')
            
            figure,
            subplot(3,1,1), plot(log.flow.t, outputs.LIPsave.p), title('com p')
            subplot(3,1,2), plot(log.flow.t, outputs.LIPsave.v, 'r'), title('com v');
            subplot(3,1,3), plot(outputs.LIPsave.p,outputs.LIPsave.v), title('phase')

        end 
        
        function plotLIP(obj, LIP) 
            log = obj.log;
            figure,
            subplot(2,2,1), plot(log.flow.t,LIP.HLIP.x), title('global pos')
            subplot(2,2,2), plot(log.flow.t,LIP.HLIP.x2f, 'r'), title('local pos');
            subplot(2,2,3), plot(log.flow.t,LIP.HLIP.dx), title('velocity')
            subplot(2,2,4), plot(log.flow.t,LIP.HLIP.u), title('step size')
        end 
        
        function plotOutputs(obj)
            log = obj.log; 
            t = log.flow.t;
            y2des = log.outputs.y2des;
            dy2des = log.outputs.dy2des; 
            outputsAct = log.outputs.y2 + y2des;
            outputsVelAct = log.outputs.dy2 + dy2des;
            
            
            figure, hold on,
            for i = 1:4 
             h = subplot(2,2,i); hold(h, 'on'); 
            plot(t, y2des(:, i),'r')
            plot(t, cumtrapz(t, dy2des(:,i)), 'b')
            legend('pos', 'vel-integral')
            end
        end 
        
        function toeVelocity(obj) 
            q = obj.log.flow.q; 
            dq = obj.log.flow.dq;
            t = obj.log.flow.t;
            N = length(q);
            leftToeVel = zeros(3, N); 
            rightToeVel = zeros(3, N); 
            
            leftToePos = zeros(3, N); 
            rightToePos = zeros(3, N); 
            
            for i = 1:N 
                qi = q(i,:)'; 
                dqi = dq(i, :)'; 
                leftToePos(:, i) = pLeftToe(qi)'; 
                leftToeVel(:, i) = vLeftToe([qi; dqi]); 
                rightToePos(:, i) = pRightToe(qi)'; 
                rightToeVel(:, i) = vRightToe([qi; dqi]); 
            end 
            
            figure, 
            subplot(2, 2, 1) 
            plot(t, leftToePos(1,:)) 
            title('left Toe p'); 
            subplot(2, 2, 2) 
            plot(t, leftToeVel(1,:)) 
            title('left Toe v'); 
            subplot(2, 2, 3) 
            plot(t, rightToePos(1,:)) 
            title('right Toe p'); 
            subplot(2, 2, 4) 
            plot(t, rightToePos(1,:)) 
            title('right Toe v'); 
            
        end 
        
    end 
    
    methods %%%% LIP analysis 
        
        % function MPCpathComp(obj) 
        %     SLIP = obj.SLIP;
        %     log = obj.log;
        %     %% check path tracking
        %     figure,
        %     subplot(2,2,1),
        %     plot(SLIP.optStepping.Xtvec_rep, SLIP.optStepping.Xvel_vec_rep,'r', log.flow.t, log.flow.dq(:,1),'b'), legend('desired', 'actual')
        %     title('X vel')
        %     subplot(2,2,2),
        %     plot(SLIP.optStepping.Xtvec_rep, SLIP.optStepping.Xpos_vec_rep,'r', log.flow.t, log.flow.q(:,1),'b'), legend('desired', 'actual')
        %     title('X global pos')
        %     subplot(2,2,3),
        %     plot(SLIP.optStepping.Ytvec_rep, SLIP.optStepping.Yvel_vec_rep,'r', log.flow.t, log.flow.dq(:,2),'b'), legend('desired', 'actual')
        %     title('Y vel')
        %     subplot(2,2,4),
        %     plot(SLIP.optStepping.Ytvec_rep, SLIP.optStepping.Ypos_vec_rep,'r', log.flow.t, log.flow.q(:,2),'b'), legend('desired', 'actual')
        %     title('Y global pos')
        % 
        %     figure, hold on,
        %     plot(SLIP.optStepping.Xpos_vec_rep, SLIP.optStepping.Ypos_vec_rep,'r')
        %     plot(log.flow.q(:,1), log.flow.q(:,2),'b'), legend('desired', 'actual')
        %     title('path')
        %     %% compare LIP and Cassie
        %     figure,
        %     h = subplot(2,4,1); hold(h, 'on'),
        %     plot(SLIP.optStepping.XFlist(2:end,1), 'r')
        %     plot(SLIP.optStepping.XLIP3vec(1,1:end-1), 'b')
        %     title('X global pos')
        %     h = subplot(2,4,2); hold(h, 'on'),
        %     plot(SLIP.optStepping.XFlist(2:end,2), 'r')
        %     plot(SLIP.optStepping.XLIP3vec(2,1:end-1), 'b')
        %     title('X local pos')
        %     h = subplot(2,4,3); hold(h, 'on'),
        %     plot(SLIP.optStepping.XFlist(2:end,3), 'r')
        %     plot(SLIP.optStepping.XLIP3vec(3,1:end-1), 'b')
        %     title('X vel')
        %     legend('act', 'LIP')
        % 
        %     h1 = subplot(2,4,4); hold(h1, 'on')
        %     plot(SLIP.optStepping.actualStepLvec(2:end), 'r')
        %     plot(SLIP.optStepping.stepLvecList(1:end-1, 2), 'b')
        %     legend('actual', 'LIP')
        %     title('step L')
        % 
        %     h = subplot(2,4,5); hold(h, 'on'),
        %     plot(SLIP.optStepping.YFlist(2:end,1), 'r')
        %     plot(SLIP.optStepping.YLIP3vec(1,1:end-1), 'b')
        %     title('Y global pos')
        %     h = subplot(2,4,6); hold(h, 'on'),
        %     plot(SLIP.optStepping.YFlist(2:end,2), 'r')
        %     plot(SLIP.optStepping.YLIP3vec(2, 1:end-1), 'b')
        %     title('Y local pos')
        %     h = subplot(2,4, 7); hold(h, 'on'),
        %     plot(SLIP.optStepping.YFlist(2:end,3), 'r')
        %     plot(SLIP.optStepping.YLIP3vec(3, 1:end-1), 'b')
        %     title('Y vel')
        %     h2 = subplot(2,4,8); hold(h2, 'on')
        %     plot(SLIP.optStepping.actualStepWvec(2:end), 'r')
        %     plot(SLIP.optStepping.stepWvecList(1:end-1, 2), 'b')
        %     legend('actual', 'LIP')
        %     title('step W')
        % end 
        
        % function MPCmodelComp(obj)
        %     %%%%%% compare the states of Cassie and the LIP in 3D projected
        %     %%%%%% in Cassie's sagittal and lateral plane. 
        %     SLIP = obj.SLIP;
        % 
        %     figure, 
        %     h = subplot(2,4,1); hold(h, 'on'),
        %     plot(SLIP.optStepping.tVec3D(2:end), SLIP.optStepping.XFlistSagittal(1, 2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.XLIP3vecSagittal(1,1:end-1), 'b')
        %     title('X global pos')
        %     h = subplot(2,4,2); hold(h, 'on'),
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.XFlistSagittal(2, 2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.XLIP3vecSagittal(2,1:end-1), 'b')
        %     title('X local pos')
        %     h = subplot(2,4,3); hold(h, 'on'),
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.XFlistSagittal(3, 2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.XLIP3vecSagittal(3,1:end-1), 'b')
        %     title('X vel')
        %     legend('act', 'LIP')
        % 
        %     h1 = subplot(2,4,4); hold(h1, 'on')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.actualStepSagittalVec(2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.stepSagittalList(1:end-1), 'b')
        %     legend('actual', 'LIP')
        %     title('step L')
        % 
        %     h = subplot(2,4,5); hold(h, 'on'),
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.YFlistLateral(1, 2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.YLIP3vecLateral(1,1:end-1), 'b')
        %     title('Y global pos')
        %     h = subplot(2,4,6); hold(h, 'on'),
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.YFlistLateral(2, 2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.YLIP3vecLateral(2, 1:end-1), 'b')
        %     title('Y local pos')
        %     h = subplot(2,4, 7); hold(h, 'on'),
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.YFlistLateral(3, 2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end),SLIP.optStepping.YLIP3vecLateral(3, 1:end-1), 'b')
        %     title('Y vel')
        %     h2 = subplot(2,4,8); hold(h2, 'on')
        %     plot(SLIP.optStepping.tVec3D(2:end), SLIP.optStepping.actualStepLateralVec(2:end), 'r')
        %     plot(SLIP.optStepping.tVec3D(2:end), SLIP.optStepping.stepLateralList(1:end-1), 'b')
        %     legend('actual', 'LIP')
        %     title('step W')
        % end 
        % 
        function checkLIPlocalApprox(obj) 
            SLIP = obj.SLIP; 
            
            %% check the LIP approximation
            figure,
            h = subplot(4,1,1); hold(h, 'on'),
            plot(SLIP.optStepping.YFlist(2:end,1), 'r')
            plot(SLIP.optStepping.Y0list(1:end-1,2), 'b')
            title('Y pos')
            h = subplot(4,1,2); hold(h, 'on'),
            plot(SLIP.optStepping.YFlist(2:end,2), 'r')
            plot(SLIP.optStepping.Y0list(1:end-1,3), 'b')
            title('Y vel')
            h = subplot(4,1,3); hold(h, 'on'),
            plot(SLIP.optStepping.XFlist(2:end,1), 'r')
            plot(SLIP.optStepping.X0list(1:end-1,2), 'b')
            title('X pos')
            h = subplot(4,1,4); hold(h, 'on'),
            plot(SLIP.optStepping.XFlist(2:end,2), 'r')
            plot(SLIP.optStepping.X0list(1:end-1,3), 'b')
            title('X vel')
            legend('act', 'est from LIP')
%             %% regress for the optimal x pos.
%             Xkp1 = SLIP.optStepping.XFlist(2:end,:);
%             Xk = SLIP.optStepping.XFlist(1:end-1, :);
%             deltaVec = zeros(2, length(Xkp1));
%             lsValue = [];
%             lsMat = [];
%             for i = 1:length(Xkp1)
%                 AminusI = SLIP.optStepping.Ahat0 - eye(2);
%                 diff = Xkp1(i,:)' - SLIP.optStepping.Ahat0*Xk(i,:)' -...
%                     SLIP.optStepping.Bhat0*SLIP.optStepping.actualStepLvec(i);
%                 deltaVec(:,i) = inv(AminusI) *diff;
%                 lsMat = [lsMat; AminusI(1,1); AminusI(2,1)];
%                 lsValue = [lsValue; diff];
%             end
%             deltaP = pinv(lsMat)*lsValue;
%             diffVec = [];
%             for i = 1:length(Xkp1)
%                 diff = [Xkp1(i,1) - deltaP; Xkp1(2)] - SLIP.optStepping.Ahat0*[Xk(i,1) - deltaP; Xk(i,2)] -...
%                     SLIP.optStepping.Bhat0*SLIP.optStepping.actualStepLvec(i);
%                 diffVec = [diffVec, diff];
%             end
        end 
        
         function checkDeadbeatLIP(obj)
           SLIP = obj.SLIP;
            log = obj.log;
            %% LIP plots  forward LIP stabilization
            set(0, 'DefaultLineLineWidth', 1);
            figure(23525), %  base plot time based
            nrow = 2; ncol = 3;
            subplot(nrow,ncol,1), plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.pelvisPos, '.r', SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.COMpos, '.b',...
                SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.HipPos, '.g')
            legend('pelvis','COM', 'hip');
            title('X position')
            subplot(nrow,ncol,2), plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.pelvisVel, '.r', SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.COMvel, '.b');
            legend('pelvis','COM');
            title('X velocity')
            
            h = subplot(nrow, ncol,3); hold(h, 'on'),
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.stepLength, '.r')
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.maxstepLength, '.b'),
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.stepLengthPreFB, '.g'),legend('stepL', 'maxStepL', 'prefeedback');
            title('X stepLength')
            
            subplot(nrow,ncol,4), plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.pelvisPos, '.r', SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.COMpos, '.b')
            legend('pelvis','COM'), title('Y position')
            subplot(nrow,ncol,5), plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.pelvisVel, '.r', SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.COMvel, '.b');
            legend('pelvis','COM'), title('Y velocity')
            
            h = subplot(nrow, ncol,6); hold(h, 'on'),
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.stepLength, '.r')
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.maxstepLength, '.b'),
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.stepLengthPreFB, '.g'),legend('stepL', 'maxStepL', 'prefeedback');
            title('Y stepLength')
            
            %%%% phase plot
            figure, hold on, box on,
            LIPx = linspace(-0.2, 0.2, 100);
            lambda =  SLIP.lambda ; %sqrt(SLIP.g/ 1) ; %
            LIP_k = SLIP.sigma_Ep;
            
            h = subplot(1, 2, 1); hold(h, 'on'), box on
            plot(SLIP.LIPstates.SSP.x.pelvisPos - SLIP.hipRoll2pelvisX/2, SLIP.LIPstates.SSP.x.pelvisVel, '-.r'),
            %LIPx = linspace(-0.2, 0.2, 100); % the quad
            %plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r'),
            xlabel('$x$', 'Interpreter','latex'),
            ylabel('$\dot{x}$', 'Interpreter','latex');
            plot(LIPx, LIP_k*LIPx, 'b',LIPx, -LIP_k*LIPx, 'b'),
            title('x phase')
            LIP_k = SLIP.sigma_En;
            
            LIPst = LIP_stabilization(SLIP.desiredVelocity(1), SLIP.EpOrEn{1});
            LIP = LIPst.setDuration(SLIP.Tssp, SLIP.Tdsp);
            %SLIP.X1Fdes = SLIP.dEn/SLIP.lambda*sinh(SLIP.Tssp*SLIP.lambda)/2;
            LIP.plotOrbit(SLIP.X1Fdes);
            
            h2 = subplot(1, 2, 2); hold(h2, 'on'), box on
            plot(SLIP.LIPstates.SSP.y.pelvisPos, SLIP.LIPstates.SSP.y.pelvisVel, '-.r'),
            %LIPx = linspace(-0.2, 0.2, 100); % the quad
            %plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r'),
            xlabel('$y$', 'Interpreter','latex'),
            ylabel('$\dot{y}$', 'Interpreter','latex');
            %plot(LIPx, LIP_k*LIPx, 'b',LIPx, -LIP_k*LIPx, 'b'),
            title('y phase')
            
            LIPst = LIP_stabilization(SLIP.desiredVelocity(2), SLIP.EpOrEn{2});
            LIP = LIPst.setDuration(SLIP.Tssp, SLIP.Tdsp);
            SLIP.X1Fdes = SLIP.dEn/SLIP.lambda*sinh(SLIP.Tssp*SLIP.lambda)/2;
            LIP.plotOrbit(SLIP.Y1Fdes);
            %%%  convergenecy
            figure,
            h = subplot(2, 2, 1); hold(h, 'on'),
            plot(SLIP.LIPstates.SSP.T, 'r')
            plot(SLIP.LIPstates.DSP.T, 'b'),title('SSP, DSP, duration')
            h = subplot(2, 2, 2); hold(h, 'on'),
            plot(SLIP.LIPstates.StepLengthX,'o')
            plot(SLIP.LIPstates.StepLengthY,'*'),
            legend('X step Length', 'Y step Length');title('actual step length')
            % feedback terms ]
            subplot(2, 2, 3),
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.Term1, '.r')
            xlabel('Time (s)'), title('feedback terms X')
            subplot(2, 2, 4),
            plot(SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.Term1, '.r')
            xlabel('Time (s)'), title('feedback terms Y')
            
            %%% understand DSP
            figure, hold on
            h = subplot(3,2,1); hold(h, 'on');
            plot(SLIP.LIPstates.SSP.x.SSPtravel, '.r'); plot(SLIP.LIPstates.DSP.x.DSPtravel, '.b') ;
            legend('SSP', 'DSP'), title('travel X dist')
            subplot(3,2,2),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.pelvisPos, '.r', ...
                SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.COMpos, '.b', ...
                SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.HipPos, '.g');
            legend('pelvis', 'Com', 'hip'), title('DSP X position');
            subplot(3,2,3),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.pelvisVel, '.r', SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.COMvel, '.b')
            legend('pelvis', 'Com'), title('DSP X velocity');
            
            h = subplot(3,2,4); hold(h, 'on');
            plot(SLIP.LIPstates.SSP.y.SSPtravel, '.r'); plot(SLIP.LIPstates.DSP.y.DSPtravel, '.b') ;
            legend('SSP', 'DSP'), title('travel Y dist')
            subplot(3,2,5),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.y.pelvisPos, '.r', ...
                SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.y.COMpos, '.b')
            legend('pelvis', 'Com'), title('DSP Y position');
            subplot(3,2,6),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.y.pelvisVel, '.r', SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.y.COMvel, '.b')
            legend('pelvis', 'Com'), title('DSP Y velocity');
            
            %%%% whole trajectories
            figure,
            subplot(2,2,1),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.pelvisPos, '.r',...
                SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.pelvisPos, '.b')
            legend('DSP', 'SSP'), title('Pelvis X pos')
            subplot(2,2,2),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.x.pelvisVel, '.r',...
                SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.x.pelvisVel, '.b')
            legend('DSP', 'SSP'), title('Pelvis X vel')
            
            subplot(2,2,3),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.y.pelvisPos, '.r',...
                SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.pelvisPos, '.b')
            legend('DSP', 'SSP'), title('Pelvis Y pos')
            subplot(2,2,4),
            plot(SLIP.LIPstates.DSP.t, SLIP.LIPstates.DSP.y.pelvisVel, '.r',...
                SLIP.LIPstates.SSP.t, SLIP.LIPstates.SSP.y.pelvisVel, '.b')
            legend('DSP', 'SSP'), title('Pelvis Y vel')
            
            %%
            figure(123321), hold on,  % phase
            plot(SLIP.pelvisYest.LIP_y_plusVec, SLIP.pelvisYest.LIP_dy_plusVec, 'o');
            LIPx = linspace(-0.2, 0.2, 100);
            plot(LIPx, SLIP.LIP_k*LIPx, 'b',LIPx, -SLIP.LIP_k*LIPx, 'b'), title('phase')
            plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r'), xlabel('$y$'), ylabel('$\dot{y}$');
            
            figure(3242), % DSP approximation , since y_DSP velocity doesn't change much so this won't be necessarily used
            h1 = subplot(2,1,1); hold(h1, 'on');
            plot(log.flow.t, log.flow.q(:, 2), 'r')
            plot(outputs.DSP_LIP_approx.t + t000, outputs.DSP_LIP_approx.y, 'b-o')
            title('position')
            h1 = subplot(2,1,2); hold(h1, 'on');
            plot(log.flow.t, log.flow.dq(:, 2), 'r')
            plot(outputs.DSP_LIP_approx.t + t000, outputs.DSP_LIP_approx.dy, 'b-o')
            title('velocity')
            
            figure(23525), %
            nrow = 2; ncol = 3;
            subplot(nrow,ncol,1), plot(SLIP.pelvisYest.t, SLIP.pelvisYest.pos, '.r', SLIP.pelvisYest.t, SLIP.pelvisYest.posNow, '.b',...
                SLIP.pelvisYest.t, SLIP.pelvisYest.COMposNow, '.k', SLIP.pelvisYest.t, SLIP.pelvisYest.LIP_y_plusVec, '.g')
            legend('est','current', 'COM', 'est_plus'), title('pelvis estimated position')
            subplot(nrow,ncol,2), plot(SLIP.pelvisYest.t, SLIP.pelvisYest.vel, '.r', SLIP.pelvisYest.t, SLIP.pelvisYest.velNow, '.b',...
                SLIP.pelvisYest.t, SLIP.pelvisYest.COMvelNow(1:end-1), '.k', SLIP.pelvisYest.t, SLIP.pelvisYest.LIP_dy_plusVec, '.g');
            legend('est','current', 'COM', 'est_plus'), title('pelvis estimated velocity')
            subplot(nrow,ncol,3), plot(SLIP.pelvisYest.t, SLIP.pelvisYest.t2impact, '.k'), title('Time to Impact estimation')
            h = subplot(nrow,ncol,4); hold(h, 'on'),
            plot(SLIP.pelvisYest.t, SLIP.pelvisYest.stepWidth, '.g'), title('stepWidth')
            plot(SLIP.pelvisYest.t, SLIP.pelvisYest.ynsfdes_lower, '.b')
            plot(SLIP.pelvisYest.t, SLIP.pelvisYest.ynsfdes_upper(1:length(SLIP.pelvisYest.t)), '.r'),
            plot(SLIP.pelvisYest.t, SLIP.pelvisYest.LIP_quad_bound, '.k')
            legend('calculated','lower','upper', 'LIP_bound')
            h = subplot(nrow,ncol,5); hold(h, 'on'),
            plot(SLIP.pelvisYest.posNow, SLIP.pelvisYest.velNow, '-.r'),
            % the quad
            LIPx = linspace(-0.2, 0.2, 100);
            plot(LIPx, SLIP.LIP_k*LIPx, 'b',LIPx, -SLIP.LIP_k*LIPx, 'b'), title('phase')
            plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r'), xlabel('$y$'), ylabel('$\dot{y}$');
            
            h = subplot(nrow,ncol,6); hold(h, 'on')
            plot( log.flow.q(:,2),log.flow.q(:,1), '.b'), title('X-Y view')
            plot(-SLIP.pelvisYest.stanceFootY, SLIP.pelvisYest.stanceFootX, 's');
            
            figure, hold on, box on,
            plot(SLIP.pelvisYest.COMposNow, SLIP.pelvisYest.COMvelNow(1:end-1), '-.r'),
            plot(SLIP.pelvisYest.posNow, SLIP.pelvisYest.velNow, '-.b'),
            
            LIPx = linspace(-0.2, 0.2, 100);% the quad
            T_SSP = 0.57;
            LIP_k = SLIP.lambda*sinh(T_SSP/2*SLIP.lambda)/cosh(T_SSP/2*SLIP.lambda);
            plot(LIPx, LIP_k*LIPx, 'b',LIPx, -LIP_k*LIPx, 'b'), title('phase') , xlabel('$y$','Interpreter', 'Latex'), ylabel('$\dot{y}$','Interpreter', 'Latex');
            plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r')
            %% forward feedback control results
            figure(23525), %
            nrow = 2; ncol = 3;
            subplot(nrow,ncol,1), plot(SLIP.pelvisXest.t, SLIP.pelvisXest.pos, '.r', SLIP.pelvisXest.t, SLIP.pelvisXest.posNow, '.b')
            legend('est','current', 'est_plus'), title('pelvis estimated position')
            subplot(nrow,ncol,2), plot(SLIP.pelvisXest.t, SLIP.pelvisXest.vel, '.r', SLIP.pelvisXest.t, SLIP.pelvisXest.velNow, '.b');
            legend('est','current', 'est_plus'), title('pelvis estimated velocity')
            subplot(nrow,ncol,3), plot(SLIP.pelvisXest.t, SLIP.pelvisXest.t2impact, 'k'), title('Time to Impact estimation')
            h = subplot(nrow,ncol,4); hold(h, 'on'),
            plot(SLIP.pelvisXest.t, SLIP.pelvisXest.stepLength, '.g'), title('stepLength')
            legend('calculated','lower','upper', 'LIP_bound')
            h = subplot(nrow,ncol,5); hold(h, 'on'),
            plot(SLIP.pelvisXest.posNow, SLIP.pelvisXest.velNow, '-.r'),
            % the quad
            LIPx = linspace(-0.2, 0.2, 100);
            plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r'), xlabel('$x$'), ylabel('$\dot{x}$');
            
            figure, hold on, box on,
            plot(SLIP.pelvisXest.posNow, SLIP.pelvisXest.velNow, '-.r'),
            plot(SLIP.pelvisXest.COMposNow, SLIP.pelvisXest.COMvelNow, '-.r'),
            
            % the quad
            LIPx = linspace(-0.2, 0.2, 100);
            T_SSP = SLIP.Walking.SSP.T;
            LIP_k = SLIP.lambda*sinh(T_SSP/2*SLIP.lambda)/cosh(T_SSP/2*SLIP.lambda);
            plot(LIPx, LIP_k*LIPx, 'b',LIPx, -LIP_k*LIPx, 'b'), title('phase') , xlabel('$x$','Interpreter', 'Latex'), ylabel('$\dot{x}$','Interpreter', 'Latex');
            plot(LIPx, SLIP.lambda*LIPx, 'r', LIPx, -SLIP.lambda*LIPx, 'r')
            
            figure,
            subplot(3,1,1), plot(SLIP.pelvisXest.t, SLIP.xLIPstatesSave.KP, '.') ; title('KP term')
            subplot(3,1,2), plot(SLIP.pelvisXest.t, SLIP.xLIPstatesSave.KI, '.') ; title('KI term')
            subplot(3,1,3), plot(SLIP.pelvisXest.t, SLIP.xLIPstatesSave.Dist_DSPx, '.') ; title('DSPx travled')
        end
        
    end 

end