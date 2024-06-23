// custom code for using ros2
#pragma once
#include <mujoco/mujoco.h>
#include <rclcpp/rclcpp.hpp>

namespace mujoco_helper
{
    // This header includes functions facilitating simulation
    // The robot would be relleased after a few seconds all nodes are launched, 
    // avoiding IK explosion when released too fast.
    // If you want to add external forces, this is also the place.
    double T_hold_in_air = 2.5;
    double T_partial_release = T_hold_in_air+1.0; 
    double T_release_horizontal = T_partial_release+2;
    std::vector<double> initial_pos{0, 0, 1};
    // it seems if make sense to provide these as serives
    static double stiffness = 1e5; 
    static double damping = 1e4; 
    static double stiffness_min = 0;
    static double damping_mid =  100; 
    static double damping_min =  1; 

    void sim_hold(mjModel *mujocoModel, mjData *mujocoData)
    {
        // Set stiffness/damping for body translation joints
        for (int i = 0; i < 2; ++i)
        {
            mujocoModel->jnt_stiffness[i] = stiffness;
            mujocoModel->dof_damping[i] = 1e4;
            mujocoModel->qpos_spring[i] = mujocoData->qpos[i];
        }
        mujocoModel->jnt_stiffness[2] = 1e3;
        mujocoModel->dof_damping[2] = 1e3;
        mujocoModel->qpos_spring[2] = mujocoData->qpos[2];
    }

    // Releases a held pelvis.
    void sim_release(mjModel *mujocoModel, mjData *mujocoData)
    {
        // Zero stiffness/damping for body translation joints
        for (int i = 0; i < 2; ++i)
        {
            mujocoModel->jnt_stiffness[i] = 0;
            mujocoModel->dof_damping[i] = 0;
        }
        mujocoModel->dof_damping[0] = 2;
        mujocoModel->dof_damping[1] = 1;

        // x direction force
        mujocoData->xfrc_applied[6] = 0;

        // // z direction force
        // mujocoData->xfrc_applied[8] =5;
        

        // // y direction periodic torque
        // double time = mujocoData->time;
        // mujocoData->xfrc_applied[10] = 100.0*sin(2*M_PI*5*time);
    }


    // main function
    void sim_helper_main(mjModel *mujocoModel, mjData *mujocoData)
    {
        double time = mujocoData->time;
        if (time <  T_hold_in_air)
            sim_hold(mujocoModel, mujocoData);
        else if (time < T_partial_release) // linear release z rpy
        {
            double rate = 1 - (time - T_hold_in_air) / (T_partial_release - T_hold_in_air);
            for (int i = 0; i < 3; ++i)
            {
                mujocoModel->jnt_stiffness[i] = stiffness * rate + stiffness_min;
                mujocoModel->dof_damping[i] = damping * rate +damping_mid;
                mujocoModel->qpos_spring[i] = mujocoData->qpos[i];
            }
        }
        else if (time < T_release_horizontal)
        {
            // linear release horizontal
            double rate = 1 - (time - T_partial_release) / (T_release_horizontal - T_partial_release);
            for (int i = 0; i < 2; ++i)
            {
                mujocoModel->jnt_stiffness[i] = 0;
                mujocoModel->dof_damping[i] = damping_mid * rate + damping_min;
             }
        }
        else // this is the time to change
            sim_release(mujocoModel, mujocoData);
   }

}