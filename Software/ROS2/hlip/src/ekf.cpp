#include "ekf.hpp"

using namespace Eigen;

namespace ekf{
    // ekf::ekf(const std::string &name) : Node(name)
    void EKF::EKF_init()
    {

        gravity_ << 0, 0, 9.81;
        C_gyro_.diagonal() << std::pow(0.1, 2), std::pow(0.1, 2), std::pow(0.1, 2);
        C_accel_.diagonal() << std::pow(3.5, 2), std::pow(3.5, 2), std::pow(3.5, 2);
        C_vo_.diagonal() << std::pow(0.01, 2), std::pow(0.01, 2), std::pow(0.01, 2), std::pow(0.01, 2);

        quternion_prior_ = VectorXd::Zero(4);
        quternion_prior_(0) = 1;

        quaternion_ = quternion_prior_;
        Cov_q_ = 0.001 * 0.001 * MatrixXd::Identity(4, 4);

        // imu_sub = create_subscription<sensor_msgs::msg::Imu>("hardware/imu_data",
        //                                                      10,
        //                                                      std::bind(&EKF::raw_data_callback, this, std::placeholders::_1));
        // publisher_filter_ = create_publisher<sensor_msgs::msg::Imu>("hardware/imu_data_filtered", 10);

        // // timer_ = create_wall_timer(std::chrono::microseconds(2500), std::bind(&ekf::timerCallback, this));
        // timer_ = create_wall_timer(std::chrono::microseconds(5000), std::bind(&ekf::timerCallback, this));
        // time_0_ = static_cast<double>(rclcpp::Clock().now().nanoseconds()) / 1e9;
    }

    // void EKF::timerCallback()
    // {
    //     // auto start1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();

    //     if (msg_num)
    //     {
    //         get_measurement();
    //         gyro_nonlinear_predict(quaternion_pred_, quaternion_, angular_vel_b_, Cov_q_, Cov_q_pred_);
    //         gyro_nonlinear_correct(quaternion_correct_, quaternion_pred_, accel_b_, Cov_q_pred_, Cov_q_correct_);

    //         quaternion_ = quaternion_correct_;
    //         Cov_q_ = Cov_q_correct_;

    //         discrete_time++;
    //     }

    //     sensor_msgs::msg::Imu imu_pub_msg;
    //     imu_pub_msg.orientation.w = quaternion_(0);
    //     imu_pub_msg.orientation.x = quaternion_(1);
    //     imu_pub_msg.orientation.y = quaternion_(2);
    //     imu_pub_msg.orientation.z = quaternion_(3);

    //     publisher_filter_->publish(imu_pub_msg);
    //     // auto stop1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
    //     // auto duration1 = static_cast<double>(stop1 - start1) / 1000000;
    //     // // std::cout << "Total loop hz elapsed: " << 1 / duration1 << " hz." << std::endl;
    //     // std::cout << 1 / duration1 << std::endl;
    // }

    void EKF::update_raw_data(Vector3d linear_acceleration, Vector3d angular_velocity, double time_now){
        // std::cout<<"start."<<std::endl;
        angular_vel_stack_.push_back(angular_velocity);
        accel_stack_.push_back(linear_acceleration);
        discrete_time_stack_.push_back(++discrete_time);
        imu_time_stack_.push_back(time_now);
        filter_quaternion_stack_.push_back(quaternion_);
        Cov_stack_.push_back(Cov_q_);
        // std::cout<<"pushed back."<<std::endl;
        auto idx_pre_ptr = std::upper_bound(imu_time_stack_.begin(),
                                            imu_time_stack_.end(), time_now); // find the first imu time bigger than vo_time_pre_
        // std::cout<<"idx_pre_ptr back."<< *idx_pre_ptr <<std::endl;
        int imu_sychron_idx = std::distance(imu_time_stack_.begin(), idx_pre_ptr) - 1;
        // std::cout<<"imu_sychron_idx back."<< imu_sychron_idx <<std::endl;
        int relative_discrete_time = discrete_time_stack_.back() - discrete_time_stack_[imu_sychron_idx];
        // std::cout<<"discrete_time_stack_ back."<<std::endl;
        // std::cout<< filter_quaternion_stack_[imu_sychron_idx](0)<<std::endl;
        // std::cout<< filter_quaternion_stack_[imu_sychron_idx](0) << filter_quaternion_stack_[imu_sychron_idx](1)<<
        //             filter_quaternion_stack_[imu_sychron_idx](2) << filter_quaternion_stack_[imu_sychron_idx](3)<<std::endl;
        quaternion_<<   filter_quaternion_stack_[imu_sychron_idx](0),filter_quaternion_stack_[imu_sychron_idx](1),
                        filter_quaternion_stack_[imu_sychron_idx](2),filter_quaternion_stack_[imu_sychron_idx](3);
        // std::cout<<"quaternion_ back."<<std::endl;
        Cov_q_ = Cov_stack_[imu_sychron_idx];

        // std::cout<<"whats wrong."<<std::endl;

        for (int i = 0; i < relative_discrete_time - 1; i++)
        {
            gyro_nonlinear_predict(quaternion_pred_, quaternion_, angular_vel_stack_[imu_sychron_idx + i], Cov_q_, Cov_q_pred_);
            gyro_nonlinear_correct(quaternion_correct_, quaternion_pred_, accel_stack_[imu_sychron_idx + i], Cov_q_pred_, Cov_q_correct_);

            if (i == 0)
            {
                quaternion_pred_ = quaternion_correct_;
                Cov_q_pred_ = Cov_q_correct_;
                // vo_nonlinear_correct(quaternion_correct_, quaternion_pred_, vo_reading, Cov_q_pred_, Cov_q_correct_);
            }

            quaternion_ = quaternion_correct_;
            Cov_q_ = Cov_q_correct_;
        }
    }

    Quaternionf EKF::predict(const Vector3d &angular_velocity){
        gyro_nonlinear_predict(quaternion_pred_, quaternion_, angular_velocity, Cov_q_, Cov_q_pred_);
        // std::cout << "after prediction:" << quaternion_pred_(0) <<" " << quaternion_pred_(1) <<" " 
        // << quaternion_pred_(2) <<" " << quaternion_pred_(3) <<  std::endl;
        Quaternionf quat(quaternion_pred_(0),quaternion_pred_(1),quaternion_pred_(2),quaternion_pred_(3));
        return quat;
    }

    Quaternionf EKF::correct(const Vector3d &acceleration){
        gyro_nonlinear_correct(quaternion_correct_, quaternion_pred_, acceleration, Cov_q_pred_, Cov_q_correct_);
        // std::cout << "after correction:" << quaternion_correct_(0) <<" " << quaternion_correct_(1) <<" " 
        // << quaternion_correct_(2) <<" " << quaternion_correct_(3) <<  std::endl;
        Quaternionf quat(quaternion_correct_(0),quaternion_correct_(1),quaternion_correct_(2),quaternion_correct_(3));
        // Quaternionf quat(1,0,0,0);
        return quat;
    }

    void EKF::get_measurement()
    {
        angular_vel_stack_.push_back(angular_vel_b_);
        accel_stack_.push_back(accel_b_);
        imu_time_stack_.push_back(imu_time_);
        discrete_time_stack_.push_back(discrete_time);
        filter_quaternion_stack_.push_back(quaternion_);
        Cov_stack_.push_back(Cov_q_);

        if (vo_new_ && imu_time_stack_.size() > 0) // log when new vo transformation was subscribed
        {
            auto start1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();

            vo_new_ = false; // meassage logged
            double _vo_time = vo_time_;
            // Vector4d vo_reading = vo_pose_quaternion_;
            // sychronize vo timestamp to the left first imu timestamp
            // imu time is the rclcpp::clock.now() when the imu_callback gets callled
            auto idx_pre_ptr = std::upper_bound(imu_time_stack_.begin(),
                                                imu_time_stack_.end(), _vo_time); // find the first imu time bigger than vo_time_pre_

            if (idx_pre_ptr == imu_time_stack_.begin())
            {
                std::cout << "not storing enough imu info, failed to correct" << std::endl;
                // should happens at the beginning of the MHE
                // discard the vo_meas
            }
            else // this if condition should only happens at during the initalization, so give some time for the imu to store enough info
            {
                int imu_sychron_idx = std::distance(imu_time_stack_.begin(), idx_pre_ptr) - 1;
                int relative_discrete_time = discrete_time_stack_.back() - discrete_time_stack_[imu_sychron_idx];
                quaternion_ = filter_quaternion_stack_[imu_sychron_idx];
                Cov_q_ = Cov_stack_[imu_sychron_idx];

                for (int i = 0; i < relative_discrete_time - 1; i++)
                {
                    gyro_nonlinear_predict(quaternion_pred_, quaternion_, angular_vel_stack_[imu_sychron_idx + i], Cov_q_, Cov_q_pred_);
                    gyro_nonlinear_correct(quaternion_correct_, quaternion_pred_, accel_stack_[imu_sychron_idx + i], Cov_q_pred_, Cov_q_correct_);

                    if (i == 0)
                    {
                        quaternion_pred_ = quaternion_correct_;
                        Cov_q_pred_ = Cov_q_correct_;
                        // vo_nonlinear_correct(quaternion_correct_, quaternion_pred_, vo_reading, Cov_q_pred_, Cov_q_correct_);
                    }

                    quaternion_ = quaternion_correct_;
                    Cov_q_ = Cov_q_correct_;
                }
            }

            auto stop1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch().count();
            auto duration1 = static_cast<double>(stop1 - start1) / 1000000;
            // std::cout << "Total loop hz elapsed: " << 1 / duration1 << " hz." << std::endl;
            std::cout << 1 / duration1 << std::endl;
        }
    }

    void EKF::gyro_nonlinear_predict(VectorXd &quaternion_pred,const VectorXd &quaternion,const Vector3d &gyro_reading, MatrixXd &Cov_q, MatrixXd &Cov_q_pred)
    {
        // std::cout << "quat input pred:" << quaternion(0) << quaternion(1) 
        // << quaternion(2) << quaternion(3) <<  std::endl;
        MatrixXd Ohm = MatrixXd::Zero(4, 4);
        gyro_2_Ohm(gyro_reading, Ohm);

        MatrixXd W = MatrixXd::Zero(4, 3);
        quat_2_W(quaternion, W);

        MatrixXd F = MatrixXd::Identity(4, 4) + dt_ / 2 * Ohm;
        // std::cout << "F matrix:\n" << F << std::endl;
        quaternion_pred = F * quaternion;

        Cov_q_pred = F * Cov_q * F.transpose() + W * C_gyro_ * W.transpose();
        // std::cout << "Convariance Q matrix:\n" << Cov_q_pred << std::endl;
        quat_norm(quaternion_pred);
    }

    void EKF::gyro_nonlinear_correct(VectorXd &quaternion_correct,const VectorXd &quaternion_pred,const Vector3d &accel_readings, MatrixXd &Cov_q_pred, MatrixXd &Cov_q_correct)
    {
        MatrixXd Rotation_pred = MatrixXd::Zero(3, 3);
        quat_2_Rot(quaternion_pred, Rotation_pred);

        Vector3d accel_hat = Rotation_pred.transpose() * gravity_;
        // std::cout << "Rotation_pred:\n" << Rotation_pred << std::endl;
        MatrixXd H = MatrixXd::Zero(3, 4);
        quat_2_H(quaternion_pred, H);

        double accel_relative_norm = accel_readings.norm() / gravity_.norm();
        // std::cout << "H matrix:\n" << H << std::endl;
        MatrixXd K = Cov_q_pred * H.transpose() * (H * Cov_q_pred * H.transpose() + accel_relative_norm * C_accel_).inverse();
        // std::cout << "K matrix:\n" << K << std::endl;
        quaternion_correct = quaternion_pred + K * (accel_readings - accel_hat);
        Cov_q_correct = (MatrixXd::Identity(4, 4) - K * H) * Cov_q_pred;

        quat_norm(quaternion_correct);
    }

    void EKF::gyro_2_Ohm(const Vector3d &gyro_reading, MatrixXd &Ohm)
    {
        // quaternion = [w,x,y,z];
        Ohm.block<1, 3>(0, 1) << -gyro_reading.transpose();

        Ohm.block<3, 1>(1, 0) << gyro_reading;

        Ohm(1, 2) = gyro_reading(2);
        Ohm(1, 3) = -gyro_reading(1);
        Ohm(2, 3) = gyro_reading(0);

        Ohm(2, 1) = -gyro_reading(2);
        Ohm(3, 1) = gyro_reading(1);
        Ohm(3, 2) = -gyro_reading(0);
    }

    void EKF::quat_2_W(const VectorXd &quaternion, MatrixXd &W)
    {
        // quaternion = [w,x,y,z];
        // W = dt /2 *[-x -y -z;
        //      w  -z  y;
        //      z   w -x;
        //      -y  x  w
        W(0, 0) = -quaternion(1);
        W(0, 1) = -quaternion(2);
        W(0, 2) = -quaternion(3);

        W(1, 0) = quaternion(0);
        W(1, 1) = -quaternion(3);
        W(1, 2) = quaternion(2);

        W(2, 0) = quaternion(3);
        W(2, 1) = quaternion(0);
        W(2, 2) = -quaternion(1);

        W(3, 0) = -quaternion(2);
        W(2, 1) = quaternion(1);
        W(2, 2) = quaternion(0);

        W = 0.5 * dt_ * W;
    }

    void EKF::quat_2_Rot(const VectorXd &quaternion, MatrixXd &Rotation)
    {

        Quaterniond q;
        q.w() = quaternion(0);
        q.x() = quaternion(1);
        q.y() = quaternion(2);
        q.z() = quaternion(3);
        Rotation = q.normalized().toRotationMatrix();
    }

    void EKF::quat_2_H(const VectorXd &quaternion, MatrixXd &H)
    {
        double w = quaternion(0);
        double x = quaternion(1);
        double y = quaternion(2);
        double z = quaternion(3);
        // std::cout << "quaternion to calculate H: " <<w<<" "<<x<<" "<<y<<" "<<z<< std::endl;
        // std::cout << "gravity to calculate H: " <<gravity_(2)<< std::endl;
        H(0, 0) = gravity_(0) * w + gravity_(1) * z - gravity_(2) * y;
        H(0, 1) = gravity_(0) * x + gravity_(1) * y + gravity_(2) * z;
        H(0, 2) = -gravity_(0) * y + gravity_(1) * x - gravity_(2) * w;
        H(0, 3) = -gravity_(0) * z + gravity_(1) * w + gravity_(2) * x;

        H(1, 0) = -gravity_(0) * z + gravity_(1) * w + gravity_(2) * x;
        H(1, 1) = gravity_(0) * y - gravity_(1) * x + gravity_(2) * w;
        H(1, 2) = gravity_(0) * x + gravity_(1) * y + gravity_(2) * z;
        H(1, 3) = -gravity_(0) * w - gravity_(1) * z + gravity_(2) * y;

        H(2, 0) = gravity_(0) * y - gravity_(1) * x + gravity_(2) * w;
        H(2, 1) = gravity_(0) * z - gravity_(1) * w - gravity_(2) * x;
        H(2, 2) = gravity_(0) * w + gravity_(1) * z - gravity_(2) * y;
        H(2, 3) = gravity_(0) * x + gravity_(1) * y + gravity_(2) * z;
        // std::cout << "H matrix:\n" << H << std::endl;
        H = 2 * H;
    }

    void EKF::quat_norm(VectorXd &quaternion)
    {
        double norm = quaternion.norm();
        quaternion = quaternion / norm;
    }

    // void EKF::raw_data_callback(const sensor_msgs::msg::Imu::SharedPtr msg)
    // {

    //     // imu_time_ = static_cast<double>(rclcpp::Clock().now().nanoseconds()) / 1e9 - time_0_; // Correctly obtaining the timestamp in seconds as a double

    //     accel_b_(0) = msg->linear_acceleration.x;
    //     accel_b_(1) = msg->linear_acceleration.y;
    //     accel_b_(2) = msg->linear_acceleration.z;

    //     angular_vel_b_(0) = msg->angular_velocity.x;
    //     angular_vel_b_(1) = msg->angular_velocity.y;
    //     angular_vel_b_(2) = msg->angular_velocity.z;
    //     msg_num++;
    // }

}