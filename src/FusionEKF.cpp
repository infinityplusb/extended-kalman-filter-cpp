#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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
    * Set the process and measurement noises
    */
    // initialise state vector
    ekf_.x_ = VectorXd(4);

    // initialise state covariance Matrix
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

    // initialise transition Matrix
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<  1, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1;

    // initialise transition Matrix
    ekf_.Q_ = MatrixXd(4, 4);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      // and initialize state.
      const float rho = measurement_pack.raw_measurements_[0];
      const float phi = measurement_pack.raw_measurements_[1];
      ekf_.x_ << rho * cos(phi), rho * sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "First measurement complete" << endl;
    return;
  }

    /**
      * Prediction
    */

    /**
      * Update the state transition matrix F according to the new elapsed time.
      * Time is measured in seconds.
      * Update the process noise covariance matrix.
      * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
    //compute the time elapsed between the current and previous measurements
    const float delta_1 = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    const float delta_2 = delta_1 * delta_1;
    const float delta_3 = delta_2 * delta_1;
    const float delta_4 = delta_3 * delta_1;

    //Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2) = delta_1;
    ekf_.F_(1, 3) = delta_1;
    // cout << "Step 01: " << endl;
    //set the acceleration noise components
    const float noise_ax = 9;
    const float noise_ay = 9;
    // cout << "Step 02: " << endl;
    //set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  delta_4/4*noise_ax, 0, delta_3/2*noise_ax, 0,
                0, delta_4/4*noise_ay, 0, delta_3/2*noise_ay,
                delta_3/2*noise_ax, 0, delta_2*noise_ax, 0,
                0, delta_3/2*noise_ay, 0, delta_2*noise_ay;
    // cout << "Step 03: " << endl;

    ekf_.Predict();
    // cout << "Predict Done!" << endl;

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // cout << "Started Radar: " << endl;
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    // cout << "Calculated the Jacobian: " << endl;
    if(!Hj_.isZero())
    {
      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    // cout << "Finished Radar: " << endl;
  }
  else {
    // cout << "Lidar: " << endl;
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
