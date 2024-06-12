#ifndef ROBOT_CONFIG_HPP
#define ROBOT_CONFIG_HPP

// #pragma once

#include <vector>
#include <string>
#include <map>

#include <Eigen/Core>
// parameters for the robot BUCKY
// things that belongs to the robot won't change

using namespace Eigen;
namespace WELL_ROBOTICS
{

    const double pi = 3.1415926;

    namespace BUCKY
    {
        static const int N_actuator_per_leg = 6;
        const int N_actuator_per_arm = 4;

        static const size_t N_motors = 20;
        static const int N_dof_constrained = 26;            // 20 + 6 floating base
        static const unsigned int N_dof_full = 30;          // 20 + 2x2 (2ankle)  + 6 floating base
        static const unsigned int N_internal_holonomic = 4; /// 4 rod distance constraints on shins
        // this would be a simplifed model
        // the actual robot has rod internal dofs that are passive.
        // this comes from the ball joints on the rod, this is not useful to model.

        // name convection better to be consistent with the robot.xml file
        // is it better to initalize these directly from the robot.xml file???
        static const std::vector<std::string> motorJointNames{
            "right_hip_roll_joint",
            "right_hip_yaw_joint",
            "right_hip_pitch_joint",
            "right_knee_joint",
            "right_inner_crank_joint",
            "right_outer_crank_joint",
            "left_hip_roll_joint",
            "left_hip_yaw_joint",
            "left_hip_pitch_joint",
            "left_knee_joint",
            "left_inner_crank_joint",
            "left_outer_crank_joint",
            "right_shoulder_yaw_joint",
            "right_shoulder_pitch_joint",
            "right_shoulder_roll_joint",
            "right_elbow_pitch_joint",
            "left_shoulder_yaw_joint",
            "left_shoulder_pitch_joint",
            "left_shoulder_roll_joint",
            "left_elbow_pitch_joint"};

        static const std::vector<std::string> rightLeg_MotorJointNames{
            "right_hip_roll_joint",
            "right_hip_yaw_joint",
            "right_hip_pitch_joint",
            "right_knee_joint",
            "right_inner_crank_joint",
            "right_outer_crank_joint"};

        static const std::vector<std::string> leftLeg_MotorJointNames{
            "left_hip_roll_joint",
            "left_hip_yaw_joint",
            "left_hip_pitch_joint",
            "left_knee_joint",
            "left_inner_crank_joint",
            "left_outer_crank_joint"};

        static const std::vector<std::string> rightArm_MotorJointNames{
            "right_shoulder_yaw_joint",
            "right_shoulder_pitch_joint",
            "right_shoulder_roll_joint",
            "right_elbow_pitch_joint"};

        const VectorXi MotorIdx_right_leg_IN_MOTORS = []()
        {
            VectorXi temp = VectorXi::LinSpaced(Sequential, 6, 0, 5); // 0, 1, --, 5
            return temp;
        }();

        const VectorXi MotorIdx_left_leg_IN_MOTORS = []()
        {
            VectorXi temp = VectorXi::LinSpaced(Sequential, 6, 6, 11); // 6, 7, --, 11
            return temp;
        }();

        const VectorXi MotorIdx_left_arm_IN_MOTORS = []()
        {
            VectorXi temp(4);
            temp << 16, 17, 18, 19;
            return temp;
        }();

        const VectorXi MotorIdx_right_arm_IN_MOTORS = []()
        {
            VectorXi temp(4);
            temp << 12, 13, 14, 15;
            return temp;
        }();

        static const std::vector<std::string> leftArm_MotorJointNames{
            "left_shoulder_yaw_joint",
            "left_shoulder_pitch_joint",
            "left_shoulder_roll_joint",
            "left_elbow_pitch_joint"};

        static const std::vector<std::string> ankleJointNames{
            "right_ankle_pitch_joint",
            "right_ankle_roll_joint",
            "left_ankle_pitch_joint",
            "left_ankle_roll_joint"};

        // this mode is simpler than the simulation model

        typedef enum
        {
            BasePosX = 0,
            BasePosY = 1,
            BasePosZ = 2,
            BaseRotX = 3,
            BaseRotY = 4,
            BaseRotZ = 5,
            // 11 joints per leg with holonomic constraints on the push rods
            // left leg  // rbdl reader reads joint order based on string names and chain .
            LeftHipRoll = 6,
            LeftHipYaw = 7,
            LeftHipPitch = 8,
            LeftInnerCrank = 9,
            LeftKnee = 10,
            LeftAnklePitch = 11,
            LeftAnkleRoll = 12,
            LeftOuterCrank = 13,

            // left arm
            LeftShoulderYaw = 14,
            LeftShoulderPitch = 15,
            LeftShoulderRoll = 16,
            LeftElbowPitch = 17,

            // right leg

            RightHipRoll = 18,
            RightHipYaw = 19,
            RightHipPitch = 20,
            RightInnerCrank = 21,
            RightKnee = 22, // first 6 are actuators
            RightAnklePitch = 23,
            RightAnkleRoll = 24,
            RightOuterCrank = 25,
            // right arm
            RightShoulderYaw = 26,
            RightShoulderPitch = 27,
            RightShoulderRoll = 28,
            RightElbowPitch = 29,
            BaseQuat = 30,
        } DofIdx; // RBDL Convention

        typedef enum
        {
            LEFT_LEG = 0,
            LEFT_ARM = 1,
            RIGHT_LEG = 2,
            RIGHT_ARM = 3
        } LIMB;
        //////////////////////////////// kinematic body points //////////////////////////
        //  get these from solidworks
        static Eigen::Vector3d pelvis_center(0, 0, 0);
        // ball joint on the crank to the origin of the crank
        static double crank_dist1 = 0.07;
        static double crank_dist2 = 0.0165;
        static Eigen::Vector3d ballJoint2innerCrank_Left(-crank_dist1, -crank_dist2, 0);
        static Eigen::Vector3d ballJoint2outerCrank_Left(-crank_dist1, crank_dist2, 0);
        static Eigen::Vector3d ballJoint2innerCrank_Right(-crank_dist1, crank_dist2, 0);
        static Eigen::Vector3d ballJoint2outerCrank_Right(-crank_dist1, -crank_dist2, 0);
        // ball joint on the foot-rotation-connector to the origin of the foot
        static Eigen::Vector3d innerBallJoint2FootJoint_Left(-0.05116, -0.0155, -0.01793);
        static Eigen::Vector3d outerBallJoint2FootJoint_Left(-0.05116, 0.0155, -0.01793);
        static Eigen::Vector3d innerBallJoint2FootJoint_Right(-0.05116, 0.0155, -0.01793);
        static Eigen::Vector3d outerBallJoint2FootJoint_Right(-0.05116, -0.0155, -0.01793);
        // foot contact points

        static double toe_dist_x = 0.0885,
                      heel_dist_x = 0.0585,
                      half_foot_width = 0.0295,
                      toe_heel_depth = 0.05193;

        static Eigen::Vector3d contactLeftInnerToe(toe_dist_x, -half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactLeftOuterToe(toe_dist_x, half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactLeftInnerHeel(-heel_dist_x, -half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactLeftOuterHeel(-heel_dist_x, half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactRightInnerToe(toe_dist_x, half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactRightOuterToe(toe_dist_x, -half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactRightInnerHeel(-heel_dist_x, half_foot_width, -toe_heel_depth);
        static Eigen::Vector3d contactRightOuterHeel(-heel_dist_x, -half_foot_width, -toe_heel_depth);

        static double shin_rod_length = 0.301;

        static Eigen::Vector3d FootPlaneCenter2AnkleCenter(0.015, 0, 0.054); // from Ankle Center to Foot Plane center

        static Eigen::Vector3d Hand2Elbow_Right(0, -0.0185, -0.24); // from elbow to hand
        static Eigen::Vector3d Hand2Elbow_Left(0, 0.0185, -0.24);   // from elbow to hand

        ///////////////////////////// joint configurations ////////////////////////////////////////////////
        // it seems cpp cannot convert string to argument names
        static std::map<std::string, size_t> motorJoint2idx{
            {"right_hip_roll_joint", RightHipRoll},
            {"right_hip_yaw_joint", RightHipYaw},
            {"right_hip_pitch_joint", RightHipPitch},
            {"right_knee_joint", RightKnee},
            {"right_inner_crank_joint", RightInnerCrank},
            {"right_outer_crank_joint", RightOuterCrank},
            {"left_hip_roll_joint", LeftHipRoll},
            {"left_hip_yaw_joint", LeftHipYaw},
            {"left_hip_pitch_joint", LeftHipPitch},
            {"left_knee_joint", LeftKnee},
            {"left_inner_crank_joint", LeftInnerCrank},
            {"left_outer_crank_joint", LeftOuterCrank},
            {"right_shoulder_yaw_joint", RightShoulderYaw},
            {"right_shoulder_pitch_joint", RightShoulderPitch},
            {"right_shoulder_roll_joint", RightShoulderRoll},
            {"right_elbow_pitch_joint", RightElbowPitch},
            {"left_shoulder_yaw_joint", LeftShoulderYaw},
            {"left_shoulder_pitch_joint", LeftShoulderPitch},
            {"left_shoulder_roll_joint", LeftShoulderRoll},
            {"left_elbow_pitch_joint", LeftElbowPitch}};

        static std::map<unsigned, unsigned int> motorJoint_idx_in_motors{
            {RightHipRoll, 0},
            {RightHipYaw, 1},
            {RightHipPitch, 2},
            {RightKnee, 3},
            {RightInnerCrank, 4},
            {RightOuterCrank, 5},
            {LeftHipRoll, 6},
            {LeftHipYaw, 7},
            {LeftHipPitch, 8},
            {LeftKnee, 9},
            {LeftInnerCrank, 10},
            {LeftOuterCrank, 11},
            {RightShoulderYaw, 12},
            {RightShoulderPitch, 13},
            {RightShoulderRoll, 14},
            {RightElbowPitch, 15},
            {LeftShoulderYaw, 16},
            {LeftShoulderPitch, 17},
            {LeftShoulderRoll, 18},
            {LeftElbowPitch, 19}};

        const VectorXi AnkleJoints_Idx_in_q = []()
        {
            VectorXi temp(4);
            temp << RightAnklePitch,
                RightAnkleRoll,
                LeftAnklePitch,
                LeftAnkleRoll;
            return temp;
        }();

        const VectorXi Right_Leg_Index_in_q = []()
        {
            VectorXi temp(8);
            temp << RightHipRoll, // 0
                RightHipYaw,      // 1
                RightHipPitch,    // 2
                RightKnee,        // 3
                RightInnerCrank,  // 4
                RightOuterCrank,  // 5
                RightAnklePitch,  // 6
                RightAnkleRoll;   // 7
            return temp;
        }();

        const VectorXi Right_Leg_main_chain_Index_in_q = []()
        {
            VectorXi temp(6);
            temp << RightHipRoll,
                RightHipYaw,     // 1
                RightHipPitch,   // 2
                RightKnee,       // 3
                RightAnklePitch, // 6
                RightAnkleRoll;  // 7
            return temp;
        }();

        const VectorXi Right_Leg_main_chain_partial_Index_in_q = []()
        {
            VectorXi temp(3);
            temp << RightHipRoll,
                RightHipPitch, // 2
                RightKnee;     // 3
            return temp;
        }();

        const VectorXi motors_idx_in_leg = []()
        {
            VectorXi out(6);
            out << 0, 1, 2, 3, 4, 5;
            return out;
        }();

        const VectorXi Left_Leg_Index_in_q = []()
        {
            VectorXi temp(8);
            temp << LeftHipRoll,
                LeftHipYaw,
                LeftHipPitch,
                LeftKnee,
                LeftInnerCrank,
                LeftOuterCrank,
                LeftAnklePitch,
                LeftAnkleRoll;
            return temp;
        }();

        const VectorXi Left_Leg_main_chain_Index_in_q = []()
        {
            VectorXi temp(6);
            temp << LeftHipRoll,
                LeftHipYaw,
                LeftHipPitch,
                LeftKnee,
                LeftAnklePitch,
                LeftAnkleRoll;
            return temp;
        }();

        // static Eigen::VectorXi motorIdx{
        // RightHipRoll,
        // RightHipYaw,
        // RightHipPitch,
        // RightKnee,
        // RightInnerCrank,
        // RightOuterCrank,
        // LeftHipRoll,
        // LeftHipYaw,
        // LeftHipPitch,
        // LeftKnee,
        // LeftInnerCrank,
        // LeftOuterCrank,
        // RightShoulderYaw,
        // RightShoulderPitch,
        // RightShoulderRoll,
        // RightElbowPitch,
        // LeftShoulderYaw,
        // LeftShoulderPitch,
        // LeftShoulderRoll,
        // LeftElbowPitch};

        static std::map<std::string, size_t> ankleJoint2idx{
            {"right_ankle_pitch_joint", RightAnklePitch},
            {"right_ankle_roll_joint", RightAnkleRoll},
            {"left_ankle_pitch_joint", LeftAnklePitch},
            {"left_ankle_roll_joint", LeftAnkleRoll}};

        struct PlaneFootContact
        { // four points
            bool right_inner_toe = false;
            bool right_outer_toe = false;
            bool right_inner_heel = false;
            bool right_outer_heel = false;

            bool left_inner_toe = false;
            bool left_outer_toe = false;
            bool left_inner_heel = false;
            bool left_outer_heel = false;
            // analyse
            bool left_whole_contact = false;
            bool right_whole_contact = false;
            unsigned int N_right_contacts = 0;
            unsigned int N_left_contacts = 0;
            // previous
            bool prev_left_whole_contact = false;
            bool prev_right_whole_contact = false;
            unsigned int prev_N_right_contacts = 0;
            unsigned int prev_N_left_contacts = 0;

            // locomotion
            bool left_touch_down = false;
            bool right_touch_down = false;
            bool left_lift_off = false;
            bool right_lift_off = false;
        };

        // PhysicalParameters
        static const double FootLength = 0.15;
        static const double FootWidth = 0.06;
        static double FootContactmomentSafety;

        static double contactFrictionCoef = 0.6; // this would be columb friction coef
        // should be a ROS parameter as well.
        static const double toeDist = 0.065;
        static const double heelDist = 0.055;

        // motor gear ratios and torque bounds,
        // may be useful, depending on the hardware
        static std::vector<double> gearRatio{
            9, 10, 9, 9, 10, 10,
            9, 10, 9, 9, 10, 10,
            6, 6, 6, 6,
            6, 6, 6, 6};

        static std::vector<double> torqueBounds{
            50, 25, 50, 50, 25, 25,
            50, 25, 50, 50, 25, 25,
            15, 15, 15, 15,
            15, 15, 15, 15};

        inline VectorXi getMotorIdx(std::vector<std::string> motorNames)
        {
            VectorXi temp = VectorXi::Zero(motorNames.size());
            for (unsigned int i = 0; i < motorNames.size(); i++)
            {
                temp(i) = BUCKY::motorJoint2idx[motorNames[i]];
            }
            return temp;
        };
    };

    typedef enum
    {
        DoNULL = 0,
        STAND = 1,
        WALK = 2
    } Locomotion_mode; // locomotion modes

    typedef enum
    {
        LEFT_STANCE = -1,
        RIGHT_STANCE = 1,
        DOUBLE_STANCE = 2,
        inAIR = 0
    } Contact_mode;

    typedef enum
    {
        x = 0,
        y = 1,
        z = 2
    } xyzidx; // xyz index

    struct WalkingCmd
    {
        double forwardVel = 0;
        double lateralVel = 0;
        double walkingPeriod = 0.35;
        double walkingHeight = 0.6;
        double stepWidth = 0.2;
        double min_stepWidth = 0.06;
        double max_stepWidth = 0.3;
        double turning_rate = 0; // rad/s yaw rate
    };
};

#endif // ROBOT_CONFIG_HPP

/*  typedef enum
      {
          BasePosX = 0,
          BasePosY = 1,
          BasePosZ = 2,
          BaseRotX = 3,
          BaseRotY = 4,
          BaseRotZ = 5,
          // 11 joints per leg with holonomic constraints on the push rods
          // right leg
          RightHipRoll = 6,
          RightHipYaw = 7,
          RightHipPitch = 8,
          RightInnerCrank = 9,
          RightOuterCrank = 10,
          RightKnee = 11, // first 6 are actuators
          RightAnklePitch = 12,
          RightAnkleRoll = 13,
          RightAnkleConnectionPitch = 14,
          RightShinInnerRodRoll = 15,
          RightShinOuterRolRoll = 16,
          // leg leg
          LeftHipRoll = 17,
          LeftHipYaw = 18,
          LeftHipPitch = 19,
          LeftInnerCrank = 20,
          LeftOuterCrank = 21,
          LeftKnee = 22,
          LeftAnklePitch = 23,
          LeftAnkleRoll = 24,
          LeftAnkleConnectionPitch = 25,
          LeftShinInnerRodRoll = 26,
          LeftShinOuterRolRoll = 27,
          // right arm
          RightShoulderYaw = 28,
          RightShoulderPitch = 29,
          RightShoulderRoll = 30,
          RightElbowPitch = 31,
          // left arm
          LeftShoulderYaw = 32,
          LeftShoulderPitch = 33,
          LeftShoulderRoll = 34,
          LeftElbowPitch = 35,
      } DofIdx_EulerXYZ; // Frost Convention


        typedef enum
      {
          BasePosX = 0,
          BasePosY = 1,
          BasePosZ = 2,
          BaseRotX = 3,
          BaseRotY = 4,
          BaseRotZ = 5,
          // 11 joints per leg with holonomic constraints on the push rods
          // right leg
          RightHipRoll = 6,
          RightHipYaw = 7,
          RightHipPitch = 8,
          RightInnerCrank = 9,
          RightOuterCrank = 10,
          RightKnee = 11, // first 6 are actuators
          RightAnklePitch = 12,
          RightAnkleRoll = 13,
          // leg leg
          LeftHipRoll = 14,
          LeftHipYaw = 15,
          LeftHipPitch = 16,
          LeftInnerCrank = 17,
          LeftOuterCrank = 18,
          LeftKnee = 19,
          LeftAnklePitch = 20,
          LeftAnkleRoll = 21,
          // right arm
          RightShoulderYaw = 22,
          RightShoulderPitch = 23,
          RightShoulderRoll = 24,
          RightElbowPitch = 25,
          // left arm
          LeftShoulderYaw = 26,
          LeftShoulderPitch = 27,
          LeftShoulderRoll = 28,
          LeftElbowPitch = 29,
          BaseQuat = 30,
      } DofIdx; // RBDL Convention

*/