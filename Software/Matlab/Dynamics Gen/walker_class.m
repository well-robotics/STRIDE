classdef walker_class < RobotLinks
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        ContactPoints
 
        mass
        nLinks

        % Center of mass position & velocity
        com_fun;
        comV_fun;
        
        LegLength_fun;
        
        LeftFoot_Rbs
        RightFoot_Rbs
        RightFootPos
        LeftFootPos
        RightAnklePos 
        LeftAnklePos
        
        LegAngle
       
        FootPos2COM  % foot position relative to COM  
        
    
        %% idx
        baseIdx = 1:3;
        motorIdx = 4:7;
        % 11 is shinspring, 12 is tarsus, 13 is heelspring
        LeftFootLinkIdx
        RightFootLinkIdx
        BaseLinkIdx = 1; 
        PelvisRot
    end
    
    methods
        function obj = walker_class(urdf)
             % Floating base model
            base = get_base_dofs('planar');
            
            % Set base DOF limits
            limits = [base.Limit];
           
            [limits.lower] = deal(-10, -10, pi/4);
            [limits.upper] = deal(10, 10, pi/4);
            [limits.velocity] = deal(20, 20, 20);
            [limits.effort] = deal(0);
            for i=1:length(base)
                base(i).Limit = limits(i);
            end
            
            obj = obj@RobotLinks(urdf,base);
       
            x = obj.States.x;
            dx = obj.States.dx;
            
            obj.com_fun = obj.getComPosition;
            obj.comV_fun = jacobian(obj.getComPosition, x) * dx;
            obj.mass = obj.getTotalMass;
             
            obj.nLinks = length(obj.Links);
            
            %% define contact frames

            r_foot_frame = obj.Joints(getJointIndices(obj, 'right_knee')); %q2_right / robot_left_knee
            obj.ContactPoints.RightToe = CoordinateFrame(...
                'Name','RightToe',...
                'Reference',r_foot_frame,...
                'Offset',[0,0,0.226],...
                'R',[pi,0,0]... % z-axis is the normal axis, so no rotation required
                );
            
            l_foot_frame = obj.Joints(getJointIndices(obj, 'left_knee')); %q2_left / robot_left_hip
            obj.ContactPoints.LeftToe = CoordinateFrame(...
                'Name','LeftToe',...
                'Reference',l_foot_frame,...
                'Offset',[0,0,0.226],...
                'R',[pi,0,0]... % z-axis is the normal axis, so no rotation required
                );
            
            %%% kneee 
            r_foot_frame = obj.Joints(getJointIndices(obj, 'right_knee'));%q2_right / robot_right_knee
            obj.ContactPoints.RightKnee = CoordinateFrame(...
                'Name','RightKnee',...
                'Reference',r_foot_frame,...
                'Offset',[0,0,0],...
                'R',[0,0,0]... % z-axis is the normal axis, so no rotation required
                );
            
            l_foot_frame = obj.Joints(getJointIndices(obj, 'left_knee'));%q2_left / robot_right_hip
            obj.ContactPoints.LeftKnee = CoordinateFrame(...
                'Name','LeftKnee',...
                'Reference',l_foot_frame,...
                'Offset',[0,0,0],...
                'R',[0,0,0]... % z-axis is the normal axis, so no rotation required
                );
            
            
            %% define other frames
            torso_frame = obj.Joints(getJointIndices(obj, 'BaseRotY'));
            obj.ContactPoints.Torso = CoordinateFrame(...
                'Name','Torso',...
                'Reference',torso_frame,...
                'Offset',[0,0,0.07],...
                'R',[0,0,0]... % z-axis is the normal axis, so no rotation required
                );
             
                 
            hipFrame = obj.Joints(getJointIndices(obj, 'BaseRotY'));
            obj.ContactPoints.HipPos = CoordinateFrame(...
                'Name','HipPos',...
                'Reference',hipFrame,...
                'Offset',[0,0,0],...
                'R',[0,0,0]... % z-axis is the normal axis, so no rotation required
                );
        end
    end

end

