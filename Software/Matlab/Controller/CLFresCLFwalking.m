classdef CLFresCLFwalking < handle
    %RESOUTPUTS 
    % calculale RES-CLF quanitities. NOTE: every single output can be
    % decoupled from the other if Q is diag (proof?)
    
    properties
        SSP
        motorIdx = 4:7;
    end
    
    methods
        function obj = CLFresCLFwalking(SSP, whichCLF)

            %% for SSP
            obj.SSP = struct('name',[],'F',[],'G',[],'Q',[],'P',[],'Pe',[],'gama',[]);
            obj.SSP.name = SSP.name;
            obj.SSP.F = [zeros(SSP.RD2size),    eye(SSP.RD2size);
                         zeros(SSP.RD2size),    zeros(SSP.RD2size)];
            obj.SSP.G = [zeros(SSP.RD2size);    eye(SSP.RD2size)];
            obj.SSP.Q = SSP.CLF_Q2;
            
            switch whichCLF
                case 'RESctleCLF' 
       
                    %%%% SSP 
                    Kp = SSP.beta^2; Kd = 2*SSP.beta;
                    A = [zeros(SSP.RD2size),eye(SSP.RD2size); -Kp*eye(SSP.RD2size), -Kd*eye(SSP.RD2size)];
                    obj.SSP.P = lyap(A',SSP.CLF_Q2);
                
                case 'REScareCLF'
                    [obj.SSP.P,   ~, ~] = care(obj.SSP.F,obj.SSP.G,obj.SSP.Q);
                   
            end
            
            % SSP RES
            e = 1/SSP.eigNum; % eigNum>1
            C3 = min(eig(obj.SSP.Q))/ max(eig(obj.SSP.P));
            emat = blkdiag(eye(SSP.RD1size),1/e.* eye(SSP.RD2size), eye(SSP.RD2size));
            obj.SSP.Pe = emat * obj.SSP.P * emat;
            obj.SSP.gama = C3/e;
            
        end
        
        function [psi0,psi1,A_hat_CLF,Lf_CLF] = CalCLF(obj, RD2, RD1, eom)
            
            % gx_hat_bottom = eom.Mass\eom.B_hat;
            % fx_bottom = eom.Mass\(eom.springVec - eom.CGvec + eom.J_hol_domain'*eom.b_Fhol2U);
            %fx_bottom = eom.Mass\(eom.springVec - eom.CGvec);
            fx_hat_bottom = eom.fx_hat_bottom;
            gx_hat_bottom = eom.gx_hat_bottom;
            
            fx = [eom.dq; fx_hat_bottom];
            
            %% set up for RD2
            %Lf2y = [dJ_y J_y]*eom.fx;  %Lf
            %LgLfy = [dJ_y J_y]*eom.gx_hat;  %A
            y2 = RD2.y;
            dy2 = RD2.dy;
            J_y2 = RD2.J_y;
            dJ_y2 = RD2.dJ_y;
            
            Lf2y = [dJ_y2, J_y2]*fx;
            LgLfy = J_y2*gx_hat_bottom;
            
            A_hat_CLF = LgLfy;
            %             Lf_CLF = Lf2y - RD2.ddy_des; % if use the derivetive term
            %             Lf_CLF = Lf2y;
            psi0 = [];
            psi1 = [];
            Lf_CLF = [];
            
            eta = [y2;dy2];
            LfV = eta'*(obj.SSP.F'*obj.SSP.Pe + obj.SSP.Pe*obj.SSP.F)*eta;
            LgV = 2*eta'*obj.SSP.Pe*obj.SSP.G;
            Veta = eta'*obj.SSP.Pe*eta;
            
            psi0 = LfV + obj.SSP.gama*Veta;
            psi1 = LgV';
            
            Lf_CLF = Lf2y;
            
            
        end
    end
end