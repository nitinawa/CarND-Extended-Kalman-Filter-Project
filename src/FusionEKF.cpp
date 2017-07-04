#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
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
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
    0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process (Q) [SN??: Seems like this is done in ProcessMeasurement below instead] and measurement noises (R) (SN: Done above)
    //SN??: Note that Q is in ekf_ (KalmanFilter class)
  */
  //SN**: initialize variables and matrices (x, F, H_laser, H_jacobian,
  //      P, etc.). Note that x, F, are members of ekf_ (KalmanFilter
  //      class) [SN??: seems like x is initilized below in
  //      ProcessMeasurement], but H_laser, H_jacobian belong to this class
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  //SN??: Hj_ << ???;  SN??: Do we need to initilize Hj? We need to update it
  //      after we get every measurement
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 0, 0, 0, 0,  //SN??: Why in 5.13, this isn't initialized
             0, 0, 0, 0,  //      with I? Since once we know the first
             0, 0, 0, 0,  //      measurement, we can assume that there
             0, 0, 0, 0;  //      is no noise.
  //ekf_.P_ << 1, 0, 0, 0,  //SN??: Why in 5.13, this isn't initialized
  //           0, 1, 0, 0,  //      with I? Since once we know the first
  //           0, 0, 1, 0,  //      measurement, we can assume that there
  //           0, 0, 0, 1;  //      is no noise.

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix (Q)
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
      
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;   //SN??: delete?

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //SN??: Why do we have to do this conversion?
      //SN: note that \rho dot is not enough to determine v_x and y_y, see Project.8 Tips and Tricks
      ekf_.x_ << measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]), measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    ekf_.Q_ = MatrixXd(4, 4);

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;  //SN??: Needed?
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
    
  ekf_.Predict(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     //SN??: Seems like we need to set H and R before calling Update()
   */

  //SN**: Seems like we should update Hj_ here
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // SN: Update H_j here
    ekf_.R_ = R_radar_;
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
