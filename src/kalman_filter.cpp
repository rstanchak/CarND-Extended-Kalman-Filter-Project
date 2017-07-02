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
    // residual
    VectorXd y = z - H_ * x_;
    // covariance
    MatrixXd S_inv = (R_ + H_ * P_ * H_.transpose()).inverse();
    // kalman gain
    MatrixXd K = P_ * H_.transpose() * S_inv;
    // state estimate
    x_ = x_ + K * y;
    // covariance estimate
    P_ = (MatrixXd::Identity(4,4) - K*H_)*P_;
}

/**
 * update the state by using Extended Kalman Filter equations
 */
void KalmanFilter::UpdateEKF(const VectorXd &z,  MeasurementFunction & h) {
    // jacobian 
    H_ = h.jacobian(x_);

    // residual
    VectorXd y = z - h.evaluate(x_);
    // covariance
    MatrixXd S_inv = (R_ + H_ * P_ * H_.transpose()).inverse();
    // kalman gain
    MatrixXd K = P_ * H_.transpose() * S_inv;
    // state estimate
    x_ = x_ + K * y;
    // covariance estimate
    P_ = (MatrixXd::Identity(4,4) - K*H_)*P_;
}
