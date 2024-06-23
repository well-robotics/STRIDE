#ifndef IK_LM_SOLVER_HPP
#define IK_LM_SOLVER_HPP

#include <iostream>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;

class IK_LM_Solver
{
public:
    IK_LM_Solver(){};
    ~IK_LM_Solver(){};
    
    void initialize(string name, unsigned int N_q, unsigned int N_y, double tolerance = 1e-6,
                    double max_iteration = 20, bool debugInfo = false)
    { // N_q is number of joints, N_y is the number of outputs.
        dim_of_q = N_q;
        dim_of_y = N_y;
        this->ik_name = name;
        
        qx_iterating.resize(N_q);
        qx_sol.resize(N_q);
        dqx_sol.resize(N_q);
        final_error.resize(N_y);
        We.resize(N_y,N_y);
        We = MatrixXd::Identity(N_y,N_y);
        // LM-methods direct
        Wn.resize(N_y,N_y);
        Wn = lambda*MatrixXd::Identity(N_y,N_y);

        this->reset(); // set them to zero

        this->error_residual_toleance = tolerance;
        this->max_iteration = max_iteration;
        this->debugInfo = debugInfo;
    };
    bool if_solved = true;

    void reset()
    { // assuming the problem is the same problem.
        qx_iterating.setZero();
        qx_sol.setZero();
        dqx_sol.setZero();
        final_error.setZero();
    };

    void solve_from_current(const VectorXd &q, const VectorXd &Ydes, const VectorXd &dYdes)
    {
        qx_iterating << q;
        this->N_iteration = IK_iteration_loop(Ydes);
        qx_sol << qx_iterating;
        error_evaluation(this->final_error, Ydes);
        MatrixXd Jacobian_of_y = MatrixXd::Zero(dim_of_y, dim_of_q);
        Jacobian_of_y = J_fun_pointer(qx_sol); 
        if (dim_of_q == dim_of_y)
            dqx_sol = Jacobian_of_y.inverse() * dYdes;
            // dqx_sol = Jacobian_of_y.completeOrthogonalDecomposition().solve(dYdes);
        else
            dqx_sol = Jacobian_of_y.completeOrthogonalDecomposition().solve(dYdes);
        // if not converged
        if (!if_solved)
        { // make sure the output is safe to execuate.
            qx_sol = q;
            dqx_sol.setZero();
        }
    }

    double error_residual_toleance = 1e-12;
    unsigned int max_iteration = 20;
    VectorXd getSolution_Position() { return qx_sol; }
    VectorXd getSolution_Velocity() { return dqx_sol; }
    // Store the solution in following variables
    VectorXd qx_sol;
    VectorXd dqx_sol;
    VectorXd final_error;
    // function pointers, output function and their Jacobian
    std::function<MatrixXd(VectorXd)> J_fun_pointer;
    std::function<VectorXd(VectorXd)> y_fun_pointer;
    //
    unsigned int N_iteration = 0;
    double computation_time = 0;
    bool debugInfo;

private:
    int dim_of_q = 1;
    int dim_of_y = 1;
    VectorXd qx_iterating; // joint angle to be solved to get to y_des
    string ik_name;
    MatrixXd We = MatrixXd::Identity(1,1);
    //For LM method
    MatrixXd Wn = MatrixXd::Identity(1,1);
    double lambda = 1e-4;
    void error_evaluation(VectorXd &error_vec, const VectorXd &Ydes)
    {
        // std::cout << ik_name <<  " act : " << y_fun_pointer(qx_iterating).transpose() << std::endl;
        // std::cout << ik_name << " des : " << Ydes.transpose() << std::endl;
        error_vec = y_fun_pointer(qx_iterating) - Ydes;
    };

    void J_error_evaluation(MatrixXd &J)
    {
        J = J_fun_pointer(qx_iterating);
    };

    unsigned int IK_iteration_loop(const VectorXd &Ydes)
    {           
        VectorXd error_function_value = VectorXd::Zero(dim_of_y);
        MatrixXd J_error = MatrixXd::Zero(dim_of_y, dim_of_q);
        double residual_error = 0.0;
        for (unsigned int i = 0; i < max_iteration; ++i)
        {
            this->error_evaluation(error_function_value, Ydes);
            residual_error = error_function_value.transpose() * error_function_value;
            if (residual_error < error_residual_toleance)
            {
                if (debugInfo)
                    std::cout << this->ik_name << " IK coverged at " << i << " iterations"
                                << "\n";
                if_solved = true;
                //std::cout<<"Jacobian\n"<<J_fun_pointer(qx_iterating)<<std::endl;
                return i;
            }

            J_error_evaluation(J_error);
            VectorXd g = J_error.transpose()*this->We*error_function_value;
            MatrixXd H = J_error.transpose()*this->We*J_error+this->Wn;
            qx_iterating = qx_iterating - (H.inverse()) * (g);
        }
        string error = this->ik_name + " IK DID NOT CONVERGE WITH ERROR" + to_string(residual_error) + "\n";
        std::cout << error;
        
        if_solved = false;
        return -1;
    }
};

#endif