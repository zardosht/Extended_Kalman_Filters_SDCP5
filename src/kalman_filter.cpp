#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  /**
   * TODO-Done: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO-Done: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO-Done: update the state by using Extended Kalman Filter equations
   */
  VectorXd y = z - radarMeasurementFunction(x_);
  // normalize the phi to be in interval (-pi, pi]
  double *phi = &y[1];
  while (*phi <= -M_PI) {
    *phi += M_PI;
  }
  while (*phi > M_PI) {
    *phi -= M_PI;
  }

  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
  
}

Eigen::VectorXd KalmanFilter::radarMeasurementFunction(const VectorXd &predictedState){
  VectorXd result = VectorXd(3);

  double px = predictedState[0];
  double py = predictedState[1];
  double vx = predictedState[2];
  double vy = predictedState[3];

  double rho = std::sqrt(std::pow(px, 2.0) + std::pow(py, 2.0));
  double phi = std::atan2(py, px);
  if (fabs(rho) < 0.0001) {
      rho = 0.0001;
   }
  double rho_dot = (px * vx + py * vy) / rho;

  result << rho, phi, rho_dot;
  return result;

}
