#include "FusionEKF.h"
#include <iostream>
#include "tools.h"
#include "Eigen/Dense"

using std::cout;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  /* Use noise_ax = 9 and noise_ay = 9 for your Q matrix. */
  var_ax_ = 9;
  var_ay_ = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // time difference
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.;

  previous_timestamp_ = measurement_pack.timestamp_;


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    VectorXd x(4);
    MatrixXd P(4, 4), F(4, 4), Q(4, 4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;
    F << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;
    Q = MatrixXd::Zero(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      // double rho_dot = measurement_pack.raw_measurements_[2];
      x(0) = rho * cos(phi);
      x(1) = rho * sin(phi);
      x(2) = 0;
      x(3) = 0;

      ekf_.Init(x, P, F, H_laser_, R_radar_, Q);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double dt4 = dt2*dt2;
  ekf_.Q_ << dt4/4*var_ax_, 0,             dt3/2*var_ax_, 0,
             0,             dt4/4*var_ay_, 0,             dt3/2*var_ay_,
             dt3/2*var_ax_, 0,             dt2*var_ax_,   0,
             0,             dt3/2*var_ay_, 0,             dt2*var_ay_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, h_radar_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
