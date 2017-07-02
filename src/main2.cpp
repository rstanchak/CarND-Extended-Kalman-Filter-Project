#include <iostream>
#include <math.h>
#include "FusionEKF.h"
#include "tools.h"

using namespace std;

int main(int argc, char* argv[]) { 
  // Create a Kalman Filter instance
  FusionEKF fusionEKF;

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  ifstream f(argv[1]);
  string sensor_measurement;

  while(getline(f, sensor_measurement))
  {
      MeasurementPackage meas_package;
      istringstream iss( sensor_measurement );
      long long timestamp;

      // reads first element from the current line
      string sensor_type;
      iss >> sensor_type;

      if (sensor_type.compare("L") == 0) {
          meas_package.sensor_type_ = MeasurementPackage::LASER;
          meas_package.raw_measurements_ = VectorXd(2);
          float px;
          float py;
          iss >> px;
          iss >> py;
          meas_package.raw_measurements_ << px, py;
          iss >> timestamp;
          meas_package.timestamp_ = timestamp;
      } else if (sensor_type.compare("R") == 0) {

          meas_package.sensor_type_ = MeasurementPackage::RADAR;
          meas_package.raw_measurements_ = VectorXd(3);
          float ro;
          float theta;
          float ro_dot;
          iss >> ro;
          iss >> theta;
          iss >> ro_dot;
          meas_package.raw_measurements_ << ro,theta, ro_dot;
          iss >> timestamp;
          meas_package.timestamp_ = timestamp;
      }
      float x_gt;
      float y_gt;
      float vx_gt;
      float vy_gt;
      iss >> x_gt;
      iss >> y_gt;
      iss >> vx_gt;
      iss >> vy_gt;
      VectorXd gt_values(4);
      gt_values(0) = x_gt;
      gt_values(1) = y_gt; 
      gt_values(2) = vx_gt;
      gt_values(3) = vy_gt;
      ground_truth.push_back(gt_values);

      //Call ProcessMeasurment(meas_package) for Kalman filter
      fusionEKF.ProcessMeasurement(meas_package);    	  

      //Push the current estimated x,y positon from the Kalman filter's state vector

      VectorXd estimate(4);

      double p_x = fusionEKF.ekf_.x_(0);
      double p_y = fusionEKF.ekf_.x_(1);
      double v1  = fusionEKF.ekf_.x_(2);
      double v2 = fusionEKF.ekf_.x_(3);

      estimate(0) = p_x;
      estimate(1) = p_y;
      estimate(2) = v1;
      estimate(3) = v2;

      estimations.push_back(estimate);

      VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);
      std::cout<<"RMSE "<<RMSE<<std::endl;
      sleep(1);


  }
}























































































