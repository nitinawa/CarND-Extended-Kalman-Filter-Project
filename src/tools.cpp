#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   TODO:
   * Calculate the RMSE here.
   */

  // SN??: Do we need rmse to have an arbitrary size according to
  //       the spatial size of estimations and ground_truth?
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0){
      cout << "Invalid estimation or ground_truth data" << endl;
      return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    
    VectorXd residual = estimations[i] - ground_truth[i];
    
    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  float rho2 = px*px+py*py;
  float rho = sqrt(rho2);
  float rho3by2 = (rho2*rho);

  //check division by zero
  if(fabs(rho2) < 0.0000001){
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/rho), (py/rho), 0, 0,
    -(py/rho2), (px/rho2), 0, 0,
    py*(vx*py - vy*px)/rho3by2, px*(px*vy - py*vx)/rho3by2, px/rho, py/rho;
  return Hj;
}
