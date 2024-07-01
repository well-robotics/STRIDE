/*
 * @class: implementation of Newton-Raphson inverse kinematics
 * @author: Yuhao Huang, Yicheng Zeng
*/

#ifndef IK_Newton_SOLVER_HPP
#define IK_Newton_SOLVER_HPP

#include <iostream>
#include <Eigen/Core>
//#include "ros_utilities/ros_utilities.hpp"

//#include "control_utilities/timer.hpp"

using namespace Eigen;
using namespace std;
// using namespace ROBOT_EOM;

class IK_Newton_Solver
{

    // IK here is defined as solving
    // delta_Y(q) = Y_act(q) - Y_desired = 0
    // delta_Y_dot(q) = Jy*dq - dY_desired = 0;
    // using newton-raphson methods to iteratively solve q
    // starting from the current configuration of q.
    // EXAMPLE A: given foot position, velocity of desired trajectory,
    //          solve leg angles q, and angular velocities dq.

    // EXAMPLE B: given foot position, velocity of desired trajectory, internal holonomic constraints
    //          solve leg angles q, and angular velocities dq.

    // assuming the problem is sqaure or N_q > Ny. what does pseudo-inverse give you?
public:
    IK_Newton_Solver(){};

    void initialize(string name,      // problem name
                    unsigned int N_q, // input dimension
                    unsigned int N_y, // output dimension
                    double tolerance = 1e-6,
                    double max_iteration = 20,
                    bool debugInfo = false)
    { // N_q is number of joints, N_y is the number of outputs.
        dim_of_q = N_q;
        dim_of_y = N_y;
        if (dim_of_q==dim_of_y){
            lambda = MatrixXd::Identity(dim_of_q,dim_of_y);
        }
        this->ik_name = name;

        qx_iterating.resize(N_q);
        qx_sol.resize(N_q);
        dqx_sol.resize(N_q);
        final_error.resize(N_y);

        this->reset(); // set them to zero
        is_initialized = true;

        this->error_residual_toleance = tolerance;
        this->max_newton_iteration = max_iteration;
        this->debugInfo = debugInfo;

        lambda = MatrixXd::Identity(N_q,N_y);
    };
    bool if_solved = true;

    ~IK_Newton_Solver(){};
    // should initialize function pointer (y_act),
    // value of and dimention of q, y

    void solve_from_current(const VectorXd &q,
                              const VectorXd &Ydes,
                              const VectorXd &dYdes)
    {
        //timer.setBegin();
        // from the current robot state
        // solve the desired q and dq
        // cout << "start from current" <<endl;
        qx_iterating << q; // start from current
        // QUESTION: should we start from the currrent configuration or previous solution
        this->N_iteration = IK_iteration_loop(Ydes);

        // if (dim_of_y == dim_of_q){
        //     this->N_iteration = IK_iteration_loop(Ydes);
        //    }
        // else
        //     IK_Opt_iteration_loop(Ydes);

        qx_sol << qx_iterating;
        
        // std::cout << "Ydes" << Ydes << std::endl; 
        // std::cout << "q_sol" << qx_sol << std::endl;
        // std::cout << qx_sol <<std::endl;
        error_evaluation(this->final_error, Ydes);
        // std::cout << this->final_error <<std::endl;
        // Compute the velocity of the qx;   J*dq = dYdes;
        MatrixXd Jacobian_of_y = MatrixXd::Zero(dim_of_y, dim_of_q);
        Jacobian_of_y = J_fun_pointer(qx_sol); // current J // q

        if (dim_of_q == dim_of_y)
            // dqx_sol = Jacobian_of_y.inverse() * dYdes;
            dqx_sol = Jacobian_of_y.completeOrthogonalDecomposition().solve(dYdes);
        else
            dqx_sol = Jacobian_of_y.completeOrthogonalDecomposition().solve(dYdes);

        // if not converged
        if (!if_solved)
        { // make sure the output is safe to execuate.
            qx_sol = q;
            dqx_sol.setZero();
        }
        // timer.setEnd(this->ik_name + " ik iteration ");
        // this->computation_time = timer.spentTime;
        // return timer.spentTime;
    };

    void reset()
    { // assuming the problem is the same problem.
        qx_iterating.setZero();
        qx_sol.setZero();
        dqx_sol.setZero();
        final_error.setZero();
        is_initialized = false;
    };

    // some constants that can be set by user
    double error_residual_toleance = 1e-12;
    unsigned int max_newton_iteration = 20;

    VectorXd getSolution_Position() { return qx_sol; }
    VectorXd getSolution_Velocity() { return dqx_sol; }
    VectorXd qx_sol;
    VectorXd dqx_sol;
    VectorXd final_error;
    // function pointers, output function and it's Jacobian
    std::function<MatrixXd(VectorXd)> J_fun_pointer;
    std::function<VectorXd(VectorXd)> y_fun_pointer;
    unsigned int N_iteration = 0;
    double computation_time = 0;
    bool debugInfo;
    //control_utilities::ExecuteTimer timer;

private:
    int dim_of_q = 1;
    int dim_of_y = 1;
    bool is_initialized = false;
    MatrixXd lambda;// = MatrixXd::Zero(1,1) ;
    VectorXd qx_iterating; // joint angle to be solved to get to y_des
    string ik_name;

    void error_evaluation(VectorXd &error_vec, const VectorXd &Ydes)
    {

        // std::cout << ik_name <<  " act : " << y_fun_pointer(qx_iterating).transpose() << std::endl;
        // std::cout << ik_name << " des : " << Ydes.transpose() << std::endl;
        error_vec = y_fun_pointer(qx_iterating) - Ydes;
        // error_vec(1) = error_vec(1)*1.5;
        // error_vec(3) = error_vec(3)*1.5;
    };

    void J_error_evaluation(MatrixXd &J)
    {
        J = J_fun_pointer(qx_iterating);
        //J(4,0) = 0;
    };

    unsigned int IK_iteration_loop(const VectorXd &Ydes)
    { // if dim_of_q = dim_of_y;
        VectorXd error_function_value = VectorXd::Zero(dim_of_y);
        MatrixXd J_error = MatrixXd::Zero(dim_of_y, dim_of_q);
        double residual_error = 0.0;
        // cout<< "start iterating, qx_init:" << qx_iterating.transpose() << endl;
        for (unsigned int i = 0; i < max_newton_iteration; ++i)
        {
            // Reinitialize the state at the current guess
            this->error_evaluation(error_function_value, Ydes);
            // Check for completion
            residual_error = error_function_value.transpose() * error_function_value; // 2 norm of the error residual
            // residual_error = error_function_value.cwiseAbs().rowwise().sum().maxCoeff(); // infinity norm
            if (residual_error < error_residual_toleance)
            {
                if (debugInfo)
                    std::cout << this->ik_name << " IK coverged at " << i << " iterations"
                              << "\n";
                if_solved = true;
                //std::cout<<"Jacobian\n"<<J_fun_pointer(qx_iterating)<<std::endl;
                return i;
            }

            // Compute the Jacobian with the new guess
            J_error_evaluation(J_error);
            // std::cout << "residual_error " << residual_error << std::endl;

            // Perform the update (Newton Raphson Method) update q with current guess of qx
            // qx_iterating = qx_iterating - J_error.completeOrthogonalDecomposition().solve(error_function_value);
            if(J_error.fullPivLu().rank() == dim_of_q){
                // std::cout << "inverse" << std::endl;
                qx_iterating = qx_iterating - (J_error.inverse()) * (error_function_value);//(J_error.inverse() - 0.1*lambda) * (error_function_value)
            }
            else{
                qx_iterating = qx_iterating - J_error.completeOrthogonalDecomposition().solve(error_function_value);
            }
            // std::cout << "new_qx: " << qx_iterating.transpose() << std::endl;
        }
        string error = this->ik_name + " IK DID NOT CONVERGE WITH ERROR" + to_string(residual_error) + "\n";
        std::cout << error;
        
        if_solved = false;
        return -1;
    };


/*
    unsigned int IK_Opt_iteration_loop(const VectorXd Ydes)
    // TODO: this does not work yet.
    { // if dim_of_q > dim_of_y; redandent joints
        // use optimization Cost = 0.5*sum(qi^2)
        // st     P(q) = Y(q) - Ydes = 0;
        // for this constrained optimization, the first order optimality condition is
        // Laguragian = cost(q) + sum(lambda_i * P_i(q))
        // J_Lagrangian_q =0,
        // J_Lagrangian_lambda = 0;
        // total equation becomes Nq + Ny;

        VectorXd error_function_value = VectorXd::Zero(dim_of_y);
        MatrixXd J_error = MatrixXd::Zero(dim_of_y, dim_of_q);
        double residual_error = 0.0;

        MatrixXd J_optimality = MatrixXd::Zero(dim_of_q + dim_of_y, dim_of_q + dim_of_y);
        VectorXd var_augmented = VectorXd::Zero(dim_of_q + dim_of_y);
        VectorXd error_augmented = VectorXd::Zero(dim_of_q + dim_of_y);

        VectorXd multiplier = VectorXd::Zero(dim_of_y);
        var_augmented << qx_iterating, // original variable to solve
            multiplier;

        for (unsigned int i = 0; i < max_newton_iteration; ++i)
        {
            // Reinitialize the state at the current guess
            this->error_evaluation(error_function_value, Ydes);
            J_error_evaluation(J_error);

            error_augmented << qx_iterating + J_error.transpose() * multiplier, // dim_of_q
                error_function_value;                                           // dim_of_y
            // Check for completion
            residual_error = error_augmented.transpose() * error_augmented; // norm of the error residual
            if (residual_error < error_residual_toleance)
            {
                std::cout << "IK coverged at " << i << " iterations"
                          << "\n";
                if_solved = true;
                return i;
            }

            // Compute the Jacobian with the new guess
            std::cout << "residual_error at " << i << " iterations " << residual_error << "\n";
            J_optimality << MatrixXd::Identity(dim_of_q, dim_of_q), J_error.transpose(),
                J_error, MatrixXd::Zero(dim_of_y, dim_of_y);

            // Perform the update (Newton Raphson Method)
            var_augmented = var_augmented - J_optimality.inverse() * error_augmented;
            qx_iterating = var_augmented.head(dim_of_q);
            multiplier = var_augmented.tail(dim_of_y);
        }
        ROS_WARN("IK DID NOT CONVERGE");
        if_solved = false;
        return -1;
    };
    */
};

#endif // IK_Newton_SOLVER_HPP
