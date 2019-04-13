#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */
#include <iostream>
using std::cout;
using std::endl;

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

void KalmanFilter::Predict() {
   /**
    * 25-8
    * KF Prediction step
    */
   x_ = F_ * x_ ;//+ u; // where u is noise
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
   /**
    * 25-8
    * KF Measurement update step
    */
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;

    // new state
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    /**
    * update the state by using Extended Kalman Filter equations
    */
    float p = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
    float p_n = (x_(0) * x_(2) + x_(1) * x_(3))/ p;
    float theta = atan2(x_(1), x_(0));

    VectorXd z_prediction = VectorXd(3);
    cout << p << endl;
    cout << theta << endl;
    cout << p_n << endl;
    z_prediction << p, theta, p_n;

    VectorXd y = z - z_prediction ;

    float angle = y(1);
    cout << M_PI << endl;
    while (angle > M_PI)
      angle -= 2* M_PI;
    while (angle < -M_PI)
      angle += 2* M_PI;
    y(1) = angle ;

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;

    //new estimate
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}
