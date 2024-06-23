

/// given the foot orientations and the calibrated joints,
// calculate the angles of the rest un-calibraed joints via Inverse kinematics
#pragma once

#include "control_common_toolbox/IK_newton.hpp"
#include "control_common_toolbox/IK_newton_sparse.hpp"
#include "control_common_toolbox/SO3.hpp"

#include "eom/robot_kinematics.hpp"
#include "control_utilities/timer.hpp"
#include <Eigen/Sparse>

#include <Eigen/Core>

#include "description/robotSensors.hpp"
using namespace Eigen;

using namespace WELL_ROBOTICS;

// should be inheriate class or friend class of anything
class IK_Calibration
{ // assuming the given kinematic functions are with full configuration as input
  // so, inside this class, the orginal kinematic functions are converged to functions
  // with only the target joint angles, that is making the problem square

    //**** unknown will be 10: 3 AK70-10 motors each leg + 2 ankle angles,
public:
    IK_Calibration(){};

    void init(ROBOT_EOM::Kinematics *kine)
    {
        kinematics = kine; // copy pointer

        this->leftCaliIdx << BUCKY::LeftHipYaw,
            BUCKY::LeftInnerCrank,
            BUCKY::LeftOuterCrank,
            BUCKY::LeftAnklePitch,
            BUCKY::LeftAnkleRoll;

        this->rightCaliIdx << BUCKY::RightHipYaw,
            BUCKY::RightInnerCrank,
            BUCKY::RightOuterCrank,
            BUCKY::RightAnklePitch,
            BUCKY::RightAnkleRoll;

        // init calibration results;
        for (int i = 0; i < leftCaliIdx.size(); ++i)
            this->q_cali.emplace(leftCaliIdx(i), 0);
        for (int i = 0; i < rightCaliIdx.size(); ++i)
            this->q_cali.emplace(rightCaliIdx(i), 0);

        init_ik_solvers(); // initialize limb ik solvers
    };

    // separate the two legs' calibrations to be robust to hardware problems; isolation
    void calibrate_from_footIMU(const VectorXd &q_robot) // ik solution that calibrate the joint angles: )
    {                                                    // get desired leg angles from the desired foot pose
        q_dummy = q_robot;

        // std::cout << "q_robot: " << q_robot << std::endl;

        // this changes the interested joints to start with a non-singular configuration.
        VectorXd qx = VectorXd::Zero(5); // all unknowns of the  ik solver

        VectorXd calib_mea(5);
        // convention:     f_output <<
        // LeftFootEulerXYZ,
        // RightFootEulerXYZ,
        // hol_rods
        // TODO add gyro foot IMU
        VectorXd calib_mea_vel = VectorXd::Zero(5);

        //*** left leg calibration
        calib_mea << kinematics->Point_LeftMidFoot.euler_angles_mes, 0, 0;

        // std::cout << "measure left euler foot: " << calib_mea.head(3).transpose() << std::endl;
        ik_leftleg_calib_solver.solve_from_current(qx,
                                                   calib_mea,
                                                   calib_mea_vel);

        //*** right leg calibration
        calib_mea << kinematics->Point_RightMidFoot.euler_angles_mes, 0, 0;

        // std::cout << "measure right euler foot: " << calib_mea.head(3).transpose() << std::endl;
        ik_rightleg_calib_solver.solve_from_current(qx,
                                                    calib_mea,
                                                    calib_mea_vel);

        // assign the solution
        for (int i = 0; i < leftCaliIdx.size(); i++)
        {
            this->q_cali[leftCaliIdx(i)] = ik_leftleg_calib_solver.qx_sol(i);
        } // assign the solution
        for (int i = 0; i < rightCaliIdx.size(); i++)
        {
            this->q_cali[rightCaliIdx(i)] = ik_rightleg_calib_solver.qx_sol(i);
        }
    }

    control_utilities::ExecuteTimer timer;

    // private:
    ROBOT_EOM::Kinematics *kinematics;
    IK_Newton_Solver ik_rightleg_calib_solver;
    IK_Newton_Solver ik_leftleg_calib_solver;

    std::map<int, double> q_cali;

    VectorXi leftCaliIdx = VectorXi::Zero(5);
    VectorXi rightCaliIdx = VectorXi::Zero(5);
    VectorXd q_dummy = VectorXd::Zero(BUCKY::N_dof_full + 1); // full configuration that is used calling the kinematics functions

private:
    void init_ik_solvers()
    {
        ik_leftleg_calib_solver.initialize("calibration left leg", 5, 5,
                                           1e-6,
                                           50);
        //**** unknown will be 10: 3 AK70-10 motors each leg + 2 ankle angles,

        // pass function pointers
        ik_leftleg_calib_solver.J_fun_pointer = [this](VectorXd q_leg)
        { return this->J_leftfoot_orientation(q_leg); };

        ik_leftleg_calib_solver.y_fun_pointer = [this](VectorXd q_leg)
        { return this->f_leftfoot_orientation(q_leg); };

        ik_rightleg_calib_solver.initialize("calibration right leg", 5, 5,
                                            1e-6,
                                            50);
        //**** unknown will be 10: 3 AK70-10 motors each leg + 2 ankle angles,

        // pass function pointers
        ik_rightleg_calib_solver.J_fun_pointer = [this](VectorXd q_leg)
        { return this->J_rightfoot_orientation(q_leg); };

        ik_rightleg_calib_solver.y_fun_pointer = [this](VectorXd q_leg)
        { return this->f_rightfoot_orientation(q_leg); };
    };

    VectorXd avoid_singularity(const VectorXd &q)
    {
        VectorXd q_non_singular = q;

        double non_singular = 0.1; // some constant away from singularity; may use a random number
        double tolerance = 0.1;
        q_non_singular(BUCKY::LeftElbowPitch) = (q(BUCKY::LeftElbowPitch) < -tolerance)
                                                    ? q(BUCKY::LeftElbowPitch)
                                                    : -non_singular; // make sure the elbow is bent
        q_non_singular(BUCKY::RightElbowPitch) = (q(BUCKY::RightElbowPitch) < -tolerance)
                                                     ? q(BUCKY::RightElbowPitch)
                                                     : -non_singular;

        q_non_singular(BUCKY::LeftKnee) = (q(BUCKY::LeftKnee) > tolerance)
                                              ? q(BUCKY::LeftKnee)
                                              : non_singular;
        q_non_singular(BUCKY::RightKnee) = (q(BUCKY::RightKnee) > tolerance)
                                               ? q(BUCKY::RightKnee)
                                               : non_singular;
        return q_non_singular;
    }

    VectorXd f_leftfoot_orientation(VectorXd &q_leg_unknown)
    {
        //**** ya = f(q)
        Vector3d lf_euler, rf_euler;
        VectorXd holonomic_rod(2);
        VectorXd f_out(5);
        //**** unknown will be 10: 3 AK70-10 motors each leg + 2 ankle angles,
        assert_size_vector(q_leg_unknown, 5);

        q_dummy(this->leftCaliIdx) = q_leg_unknown; // use the current leg angles

        lf_euler << kinematics->get_body_eulerAngles(q_dummy, kinematics->Point_LeftMidFoot);

        holonomic_rod << kinematics->getClosedLoopHolonomicEachLeg(q_dummy, BUCKY::LEFT_LEG),

            f_out << lf_euler, holonomic_rod;
        return f_out;
    };

    MatrixXd J_leftfoot_orientation(VectorXd &q_leg_unknown)
    {
        // matrix to storge Jfoot
        MatrixXd J_foot_pose_robot(6, BUCKY::N_dof_full);
        J_foot_pose_robot.setZero();

        assert_size_vector(q_leg_unknown, 5);

        q_dummy(this->leftCaliIdx) = q_leg_unknown; // swap

        MatrixXd J_lf_euler(3, 5);

        // calculate left foot pose jacobain
        kinematics->get_Jacobian_BodyPoint6d_Euler(q_dummy,
                                                   kinematics->Point_LeftMidFoot,
                                                   J_foot_pose_robot);

        J_lf_euler << J_foot_pose_robot(Eigen::indexing::seq(0, 2), this->leftCaliIdx);

        // get the push rods jacobian
        MatrixXd J_hol_left = kinematics->getClosedLoopJacobianEachLeg(q_dummy, BUCKY::LEFT_LEG);

        MatrixXd J_leg(5, 5);
        J_leg << J_lf_euler,
            J_hol_left(Eigen::indexing::all, this->leftCaliIdx);

        return J_leg;
    };

    VectorXd f_rightfoot_orientation(VectorXd &q_leg_unknown)
    {
        //**** ya = f(q)
        Vector3d rf_euler;
        VectorXd holonomic_rod(2);
        VectorXd f_out(5);
        //**** unknown will be 10: 3 AK70-10 motors each leg + 2 ankle angles,
        assert_size_vector(q_leg_unknown, 5);

        q_dummy(this->rightCaliIdx) = q_leg_unknown; // use the current leg angles

        rf_euler << kinematics->get_body_eulerAngles(q_dummy, kinematics->Point_RightMidFoot);

        holonomic_rod << kinematics->getClosedLoopHolonomicEachLeg(q_dummy, BUCKY::RIGHT_LEG);

        f_out << rf_euler, holonomic_rod;
        return f_out;
    };

    MatrixXd J_rightfoot_orientation(VectorXd &q_leg_unknown)
    {
        // matrix to storge Jfoot
        MatrixXd J_foot_pose_robot(6, BUCKY::N_dof_full);
        J_foot_pose_robot.setZero();

        assert_size_vector(q_leg_unknown, 5);

        q_dummy(this->rightCaliIdx) = q_leg_unknown; // swap

        MatrixXd J_rf_euler(3, 5);

        // calculate right foot pose jacobain
        kinematics->get_Jacobian_BodyPoint6d_Euler(q_dummy,
                                                   kinematics->Point_RightMidFoot,
                                                   J_foot_pose_robot);

        J_rf_euler << J_foot_pose_robot(Eigen::indexing::seq(0, 2), this->rightCaliIdx);

        // get the push rods jacobian
        MatrixXd J_hol_right = kinematics->getClosedLoopJacobianEachLeg(q_dummy, BUCKY::RIGHT_LEG);

        MatrixXd J_leg(5, 5);
        J_leg << J_rf_euler,
            J_hol_right(Eigen::indexing::all, this->rightCaliIdx);

        return J_leg;
    };
};
