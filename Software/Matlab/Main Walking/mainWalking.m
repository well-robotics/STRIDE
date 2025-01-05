clc, clear
%% path 
addpath(genpath('D:\Working\packages\bipedal robotics\Matlab')); % change this path to where you put the matlab folder
addpath(genpath('D:\Working\packages\calculation\qpOASES')); % change this path to where you install qpoases


%% initial configturation 
Frequency =  1e+3; %1e+3
ifubound = 1;  % 1: use torque bounds
qsize = 7; 
miu = 0.6;  % friction coeff
vxdes = 0.1; % x velocity desired
COMheight = 0.33;
tSSP = 0.4;  
swingHeight = 0.02;  
impactVel = -0.1; % 
torsoAngle = 0; %   
gaitType = 'P1';  %%%% 'P1', 'P2' 
stepL_left = -0.05;  %%%% one step size for P2 orbits 

q0 = [0; 0; torsoAngle;  
      pi/2; -pi;  
      -pi/2; pi]; 
stanceFoot = pLeftToe(q0); 
q0(2) = -stanceFoot(2); 
dq0 = zeros(7,1); 
x0 = [q0; dq0]; 
leftFootOffset = [0; 0]; %%% for walking 
q0 = InverseKinematicsWalker(COMheight, q0, torsoAngle, leftFootOffset);
% showWalker(q0); 
comPosZ = COMPosition(q0)*[0;0;1]; 
param = struct('COMz', comPosZ, 'tSSP', tSSP,  'vxdes', vxdes,  'swingHeight', swingHeight,...
                'torsoAngle', 0, 'impactVel', impactVel, 'gaitType', gaitType, 'stepL_left', stepL_left);
contact = AmberConstants.RightFootContact;
%%%%%%%
QPvariableType  =  'OSC';  %%% 'onlyU' 'U-Fhol', 'U-Fhol-ddq', 'OSC'
QP_p = 1e+5; % parameters in CLF relaxation cost

period = 1/Frequency; 
log = struct(); 
log.finishingIdx = []; 
log.flow = struct('t', [], 'q', [], 'dq',[],  'u', [], 'uDM', [], 'comX', []);
log.QP = struct('F', [], 'delta' , [], 'exitflag', [], 'numiter',[], 'fval', [], 'GRF', []);
log.DMC = struct('u', [], 'cost', [], 'F', [], 'exitflag', []);
log.outputs = struct('y2', [], 'dy2', [], 'y2des', [], 'dy2des', [], 'yNext', []);
log.param = param; 
logOpt = struct('alpha', [], 'p',[]);
outputs = GenOutputsFiveLinkWalker(logOpt, tSSP, comPosZ, swingHeight, impactVel, vxdes, gaitType, stepL_left, period); 

outputs.qPelvis0 = q0(3);
eomType = struct('QPvariableType', QPvariableType);
eom = EOM_walker(7, eomType); 
 
eventSim = customSimWalking2D(contact); 
qp = QP_Walking(miu, QPvariableType);
q = q0; 
dq = dq0; 
d = 1; 
Ncycles = 8;  % 4
Nperiods = round(tSSP/period)*Ncycles;
t0 = 0; 
lastContact = contact;
ddq = zeros(7,1); 
QPfail = false; 

%%
tic
u0 = []; 
odeopts = odeset('MaxStep', 1e-3, 'RelTol', 1e-4, 'AbsTol',1e-4);
tIdx = 0;
tIdx0 = 1;

while (d<= Ncycles && tIdx <Nperiods)
    for tIdx = tIdx0:Nperiods
        tIdx
        if eventSim.contact ~= lastContact
            outputs.SSPnsf0 = [];
            contact = eventSim.contact;
            lastContact = eventSim.contact;
            log.finishingIdx = [log.finishingIdx, tIdx];
            d = d + 1;
             tIdx0 = tIdx + 1;
            break;
        end
        if q(2) < 0 %%% falling
            break
        end
        eom.updateWalking(q, dq, contact);
        outputs = outputs.getOutputs(q, dq, ddq, tIdx - tIdx0+1, contact);
 
        [u, exitflag, F_GRF, delta, fval, numiter, QPfail] = qp.constructAndSolve(contact, u0, eom, outputs);
        [tt, x, ddq] = eventSim.sim(period, [q; dq], eom, u, odeopts, t0, QPfail);
        xnext = x(end, :);
        q  = xnext(1:qsize)'; 
        dq = xnext(qsize+1:end)';
        tt = tt + t0*ones(size(tt));
        t0 = tt(end);
        
        % data logging
        log.flow.t      =  [log.flow.t;     t0];
        log.flow.q      =  [log.flow.q;     q'];
        log.flow.dq     =  [log.flow.dq;   dq'];
        log.flow.u      =  [log.flow.u;     u'];
        log.flow.comX   =  [log.flow.comX, [outputs.comPos; outputs.comVel]]; 
        log.QP.GRF = [log.QP.GRF, F_GRF]; 
        log.outputs.yNext = [log.outputs.yNext;   outputs.y_kp1'];
        log.outputs.y2 = [log.outputs.y2;         outputs.RD2.y_unscaled'];
        log.outputs.dy2 = [log.outputs.dy2;       outputs.RD2.dy_unscaled'];
        log.outputs.y2des = [log.outputs.y2des;   outputs.RD2.y_des_unscaled'];
        log.outputs.dy2des = [log.outputs.dy2des; outputs.RD2.dy_des_unscaled'];
    end

end
toc

%%
anasim = analysisWalker(log, outputs, Ncycles, eventSim, Frequency);
anasim.plotTracking(); 
anasim.animate();
%%%%%%%% export trajectory %%%%%%%%%%%
% writematrix(anasim.log.QP.GRF','fs.txt','Delimiter','space','WriteMode','overwrite')
% writematrix(anasim.log.outputs.y2des,'ydes.txt','Delimiter','space','WriteMode','overwrite'); 
% writematrix(anasim.log.outputs.dy2des,'dydes.txt','Delimiter','space','WriteMode','overwrite'); 
% writematrix(anasim.log.flow.u,'u.txt','Delimiter','space','WriteMode','overwrite')

anasim.plotOutputs(); 















