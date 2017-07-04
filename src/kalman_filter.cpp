#include "kalman_filter.h"
//SN
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;  //SN??: How is this initialize? Its dimension is going to keep
              //      chaning from 2x4 and 3x4 for RADAR and LIDAR
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict(float delta_T) {
  /**
  TODO:
    * predict the state
  */
  float dt = delta_T;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
    
  //Modify the F matrix so that the time is integrated
  //F_ = MatrixXd(4, 4);  SN: This is done in FusionEKF()
  F_ << 1, 0, dt, 0,
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;
    
  //set the process covariance matrix Q
  //float noise_ax = 9.0*9.0;  // P_0 = zeros(4, 4), RMSE = 0.1984, 0.8885, 0.9338, 3.7263
  //float noise_ay = 9.0*9.0;  // P_0 = eyes(4, 4), RMSE = 0.1956, 0.8894, 0.9453, 3.7365
  float noise_ax = 9.0;  // P_0 = zeros(4, 4), RMSE = 0.1707, 0.6670, 0.6269, 1.6165
  float noise_ay = 9.0;  // P_0 = eyes(4, 4), RMSE = 0.1403, 0.6662, 0.5804, 1.6363
  //Q_ = MatrixXd(4, 4);  SN: This is done in
  //FusionEKF::ProcessMeasurement()
  Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
    0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
    dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
    0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
    
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  //SN1: Compute error y = z - H_*x_
  VectorXd y(2);
  y << 0, 0;
  y = z - H_*x_;
  //SN2: Compute the Kalman gain using H_ = Hj_ and R_
  Eigen::MatrixXd S;
  Eigen::MatrixXd K;
  S = MatrixXd(2, 2);
  K = MatrixXd(4, 2);
  S = H_*P_*H_.transpose() + R_;
  //SN??: Do we need to check if S is invertible?
  K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  Eigen::MatrixXd I;
  I = MatrixXd(4, 4);
  I << 1.0, 0, 0, 0,
    0, 1.0, 0, 0,
    0, 0, 1.0, 0,
    0, 0, 0, 1.0;
  P_ = (I - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  //SN1: Compute error y = z - h(x_)
  VectorXd y(3);
  y << 0,0,0;
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho = sqrt(px*px+py*py);
    
  //check division by zero
  if(fabs(rho) < 0.0000001){
      std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return;
  }
  y(0) = z(0) - rho;
  y(1) = z(1) - atan2(py, px);
  // SN**: Need to check the range of argument of atan
  while(y(1) < -3.14159265358979323846){
      y(1) = y(1) + 2*3.14159265358979323846;
  }
  y(2) = z(2) - (px*vx + py*vy)/rho;
    
  //SN2: Compute the Kalman gain using H_ = Hj_ and R_
  Eigen::MatrixXd S;
  Eigen::MatrixXd K;
  S = MatrixXd(3, 3);
  K = MatrixXd(4, 3);
  S = H_*P_*H_.transpose() + R_;
  //SN??: Do we need to check if S is invertible?
  K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  Eigen::MatrixXd I;
  I = MatrixXd(4, 4);
  I << 1.0, 0, 0, 0,
       0, 1.0, 0, 0,
       0, 0, 1.0, 0,
       0, 0, 0, 1.0;
  P_ = (I - K*H_)*P_;
}
