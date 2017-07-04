#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

/**
 * predict the state
 */
void KalmanFilter::Predict() {
    std::cout<<"<Predict>"<<std::endl;
    std::cout<<"x_:"<<x_<<std::endl;
    std::cout<<"P_:"<<P_<<std::endl;
    std::cout<<"F_:"<<F_<<std::endl;
    std::cout<<"Q_:"<<Q_<<std::endl;
   // state prediction
   x_ = F_ * x_;

   // covariance prediction
   P_ = F_ * P_ * F_.transpose()  + Q_;
    std::cout<<"x_:"<<x_<<std::endl;
    std::cout<<"P_:"<<P_<<std::endl;
    std::cout<<"</Predict>"<<std::endl;
}

/**
 * update the state by using Kalman Filter equations
 */
void KalmanFilter::Update(const VectorXd &z) {
    std::cout<<"<Update>"<<std::endl;
    std::cout<<"x_:"<<x_<<std::endl;
    std::cout<<"P_:"<<P_<<std::endl;
    std::cout<<"z:"<<z<<std::endl;
    std::cout<<"H_:"<<H_<<std::endl;
    std::cout<<"R_:"<<R_<<std::endl;
    // residual
    VectorXd y = z - H_ * x_;
    std::cout<<"y:"<<y<<std::endl;
    // covariance
    MatrixXd S_inv = (R_ + H_ * P_ * H_.transpose()).inverse();
    std::cout<<"S_inv:"<<S_inv<<std::endl;
    // kalman gain
    MatrixXd K = P_ * H_.transpose() * S_inv;
    std::cout<<"K:"<<K<<std::endl;
    // state estimate
    x_ = x_ + K * y;
    // covariance estimate
    P_ = (MatrixXd::Identity(4,4) - K*H_)*P_;
    std::cout<<"x_:"<<x_<<std::endl;
    std::cout<<"P_:"<<P_<<std::endl;
    std::cout<<"</Update>"<<std::endl;
}

/**
 * update the state by using Extended Kalman Filter equations
 */
void KalmanFilter::UpdateEKF(const VectorXd &z,  MeasurementFunction & h) {
    std::cout<<"<UpdateEKF>"<<std::endl;
    std::cout<<"x_:"<<x_<<std::endl;
    std::cout<<"P_:"<<P_<<std::endl;
    // jacobian 
    H_ = h.jacobian(x_);
    std::cout<<"z:"<<z<<std::endl;
    std::cout<<"H_:"<<H_<<std::endl;
    std::cout<<"R_:"<<R_<<std::endl;

    // residual
    VectorXd y = h.residual( z, x_ );
    std::cout<<"y:"<<y<<std::endl;
    // covariance
    MatrixXd S = (R_ + H_ * P_ * H_.transpose());
    std::cout<<"S:"<<S<<std::endl;
    // kalman gain
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    std::cout<<"K:"<<K<<std::endl;
    // state estimate
    x_ = x_ + K * y;
    // covariance estimate
    P_ = (MatrixXd::Identity(4,4) - K*H_)*P_;
    std::cout<<"x_:"<<x_<<std::endl;
    std::cout<<"P_:"<<P_<<std::endl;
    std::cout<<"</UpdateEKF>"<<std::endl;
}
