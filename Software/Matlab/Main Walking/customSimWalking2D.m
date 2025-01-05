classdef customSimWalking2D < handle

    %%%%% 
    properties (Access = public)
        %
        contact
        precontact % previous contact 
        leftFootForce

        rightFootForce
        holonomic_forces
        J_hol
        dJ_hol
        contactSequence
        contactTime
        
        guard
        countdown = 10;
        preImpact
        postImpact
        
        ndof = 7
        vho % contains vho state % a guard is needed.
        
        relabel = [eye(3), zeros(3,4); 
                        0, 0, 0, 0, 0, 1, 0; 
                        0, 0, 0, 0, 0, 0, 1; 
                        0, 0, 0, 1, 0, 0, 0; 
                        0, 0, 0, 0, 1, 0, 0]
        z_terrain
        
             
        %%%% disturbance analysis
        disturbance = struct('t', [], 'Fx', [], 'Fy', [])
        PushForce = 0; 
        Fx = [];
        Fy = []; 
        ifaddDisturb = 0;
    end

    methods
        function obj = customSimWalking2D(contact)
            obj.contact = contact;
            obj.leftFootForce = [];
            obj.rightFootForce = [];

            obj.holonomic_forces = [];
            obj.contactSequence = [obj.contact];
            obj.contactTime = 0;
            
            obj.guard = struct('time',[],'force',[],'pos', []);
            obj.countdown = 10; % hack to make sure landing is triggered after flight for a while
            obj.preImpact = [];
        end
    end
    
    methods
        
        function [T, X, ddq] = sim(obj, period, x0, eom, u, odeopts, tNow, QPfail)
            t0 = 0;
            q_size = obj.ndof; 
            % if QP fails and force is negative then step to next DSPain
            % because event function is not going to be triggered in this
            % situation
            obj.calHolonomic_forces(x0, eom, q_size, u)
            
            while abs(t0 - period) ~= 0
                options = odeset(odeopts,'Events',@(t,x) obj.guardsSym(t,x, eom, q_size, u,tNow));
                [T, X, ~, ~, ie] = ode45(@(t,x) obj.DynamicsAsym(t, x, obj.contact, eom, q_size, u), [t0,period], x0, options)
                t0 = T(end);
                x0 = X(end,:)';
                if (t0 - period) ~= 0 && obj.countdown == 0 % guard being reached 
                    switch obj.contact
                        case AmberConstants.RightFootContact % 'RightFootContact'
                            obj.contact = AmberConstants.LeftFootContact; %'LeftFootContact';
                            obj.preImpact = struct('q', x0(1:obj.ndof), 'dq', x0(obj.ndof+1:end));
                            [x0, Fimp] = obj.LandingImpact(x0, eom, q_size); % why the impact force is not correct cause it is impact impluse
                            obj.postImpact = struct('q', x0(1:obj.ndof), 'dq', x0(obj.ndof+1:end));
                            obj.calHolonomic_forces(x0, eom, q_size, u); 
                            obj.countdown = round( 1/period /10); % avoid continuous trigger of switchings
                        case AmberConstants.LeftFootContact %'LeftFootContact'
                           obj.contact = AmberConstants.RightFootContact; % 'RightFootContact';
                            obj.preImpact = struct('q', x0(1:obj.ndof), 'dq', x0(obj.ndof+1:end));
                            [x0, Fimp] = obj.LandingImpact(x0, eom, q_size);
                            obj.postImpact = struct('q', x0(1:obj.ndof), 'dq', x0(obj.ndof+1:end));
                            obj.calHolonomic_forces(x0, eom, q_size, u); 
                            obj.countdown = round( 1/period /10);   % hack
                        case  AmberConstants.DoubleSupport %'DoubleSupport'  %%%% should not go into here 
                            disp('wrong in sim')
                            if ie == 1
                                obj.contact = AmberConstants.RightFootContact; %'RightFootContact';
                            elseif ie == 2
                                obj.contact = AmberConstants.RightFootContact; %'LeftFootContact';
                            end
                            obj.countdown = round( 1/period /10);    
                    end
                    obj.contactSequence = [obj.contactSequence,  obj.contact];
                    obj.contactTime = [obj.contactTime, tNow+t0];
                end
            end
            dX =  DynamicsAsym(obj, t0, x0,  obj.contact, eom, q_size, u);
            ddq = dX(q_size+1:end); 
            if obj.countdown > 0
                obj.countdown = obj.countdown - 1;
            end

            % update contact forces
            obj.calHolonomic_forces(x0, eom, q_size, u)
            switch obj.contact
                case AmberConstants.RightFootContact % 'RightFootContact'
                    obj.leftFootForce = [obj.leftFootForce, zeros(2,1)];
                    obj.rightFootForce = [obj.rightFootForce, obj.holonomic_forces(1:2)];
                case AmberConstants.LeftFootContact  %'LeftFootContact'
                    obj.leftFootForce = [obj.leftFootForce, obj.holonomic_forces(1:2)];
                    obj.rightFootForce = [obj.rightFootForce, zeros(2,1)];
                case AmberConstants.DoubleSupport % 'DoubleSupport'
                    obj.leftFootForce = [obj.leftFootForce, obj.holonomic_forces(1:2)];
                    obj.rightFootForce = [obj.rightFootForce, obj.holonomic_forces(3:4)];
            end
        end

        function  xdot = DynamicsAsym(obj, t, x, contact, eom, q_size, u)
            q = x( 1:q_size );
            dq = x( q_size+1:2*q_size );
            eom.updateWalking(q,dq, obj.contact ); 
            obj.J_hol = eom.J_hol_domain;
            obj.dJ_hol = eom.dJ_hol_domain;
            if obj.ifaddDisturb == 1
                JpelvisPos = J_pelvis_pos(q);
                h = eom.h - transpose(JpelvisPos(1:2,:))*[obj.Fx; obj.Fy];
            else
                h = eom.h;
            end
            
            XiInv = obj.J_hol * (eom.Mass \ transpose(obj.J_hol));
            eomRightsidepart =  - h + eom.Bu*u;
            obj.holonomic_forces = -XiInv \ (obj.dJ_hol * dq + obj.J_hol * (eom.Mass \eomRightsidepart));
            
            eomRightside = eomRightsidepart + obj.J_hol'*obj.holonomic_forces;
            xdot = [dq;
                    eom.Mass\eomRightside];
        end
        
        function [value,isterminal,direction] = guardsSym(obj, t, x, eom, q_size, u, t0)
            q = x(1:q_size);
            dq = x(q_size+1:2*q_size);
           %  eom.updateWalking(q,dq, obj.contact ); 
            
            switch obj.contact
                case AmberConstants.DoubleSupport % 'DoubleSupport'
                    disp('wrong in sim'); 
                    value(1) = obj.holonomic_forces(2); % Fz on left foot
                    value(2) = obj.holonomic_forces(4); 
                    obj.guard.time = [obj.guard.time,t0+t];
                    obj.guard.force = [obj.guard.force,value];
                    obj.guard.pos = [obj.guard.pos,0];
                    isterminal(1) = 1;
                    direction(1) = -1;
                    isterminal(2) = 1;
                    direction(2) = -1;
                case AmberConstants.LeftFootContact
                    value = [0,1]*pRightToe(q);
                    obj.guard.time = [obj.guard.time,t0+t];
                    obj.guard.pos = [obj.guard.pos,value];
                    obj.guard.force = [obj.guard.force,0];
                    isterminal = 1;
                    direction = -1;
                case AmberConstants.RightFootContact
                    value = [0,1]*pLeftToe(q);
                    obj.guard.time = [obj.guard.time,t0+t];
                    obj.guard.pos = [obj.guard.pos,value];
                    obj.guard.force = [obj.guard.force,0];
                    isterminal = 1;
                    direction = -1;
            end
        end

        function [xpost, Fimp] = LandingImpact(obj,xpre,eom,q_size)
            % post impact with two feet
            % D(dqp - dqm) = J'*Fimp;
            %          J*dqp = 0
            % set up fo\dot{q}  linear system equations
            % closed form solution 
            % dq_plus = (I - inv(M)*J'*inv(J*inv(M)*J')*J) * dq_minus;
            % impF = - inv(J*inv(M)*J')*J*dq_minus
            q = xpre(1:q_size);
            dq = xpre(q_size+1:2*q_size);
            eom.updateWalking(q,dq, obj.contact ); 
            
            %%% 
            A = [eom.Mass, -transpose(eom.J_hol_domain );
                eom.J_hol_domain, zeros(size(eom.J_hol_domain, 1))];
            b = [eom.Mass*dq; zeros(size(eom.J_hol_domain,1),1)];
            sol_lin = A\b;
            dq_post = sol_lin(1:q_size);
            Fimp = sol_lin(q_size+1:end);
            
%             %%%% relabel 
%             q = obj.relabel*q; 
%             dq_post = obj.relabel*dq_post; 
            xpost = [q;dq_post];
        end
        
        function calHolonomic_forces(obj, x, eom, q_size, u)
            q = x(1:q_size);
            dq = x(q_size+1:2*q_size);
            eom.updateWalking(q, dq, obj.contact ); 
            
            obj.J_hol = eom.J_hol_domain; 
            obj.dJ_hol = eom.dJ_hol_domain;
            
            XiInv = obj.J_hol * (eom.Mass \ transpose(obj.J_hol));
            eomRightsidepart =  - eom.CGvec + eom.Bu*u;
            obj.holonomic_forces = -XiInv \ (obj.dJ_hol * dq + obj.J_hol * (eom.Mass \eomRightsidepart)); 
            
        end
        
        % function AddDisturbance(obj, tNow) %%%% %
        %     %% TODO start for a whole step, instead of a time duration
        %     %%% add disturbance
        %     if (tNow > 3.833 && tNow < 4.21 )  % ||(tNow > 14 && tNow < 14.1 ) || (tNow > 21 && tNow < 21.1 )
        %         obj.disturbance.t = [obj.disturbance.t, tNow];
        %         obj.disturbance.Fx = [obj.disturbance.Fx, obj.Fx];
        %         obj.disturbance.Fy = [obj.disturbance.Fy, obj.Fy];
        %         obj.ifaddDisturb = 1;
        %     else
        %         obj.Fx = obj.PushForce; %33*5;
        %         obj.Fy = 0; %33*5;
        %         obj.ifaddDisturb = 0;
        %    end
        % end
    end

    
end

