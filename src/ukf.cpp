#include "ukf.h"
#include "Eigen/Dense"
#include  <iostream>

using namespace std;


using Eigen::MatrixXd;
using Eigen::VectorXd;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
 
  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;

  // State Dimension
  n_x_ = 5;

  // Argument Dimension
  n_aug_ = 7;
  
  // Measurement Dimension 
  n_z_laser_ = 2;
  n_z_radar_ = 3;
  
  // spreding Parameter
  lambda_ = 3.0-n_aug_;
  //lambda_ = 3.0-n_x_;

  // Sigma Points
  Xsig_pred_ = MatrixXd(n_x_ , 2*n_aug_+1);
  //Weigths
  weights_ = VectorXd(2*n_aug_ + 1 );

  // time when the state is true, in us
  time_us_ = 0.0;
  
  // InitParameter
  x_.setZero();
  Xsig_pred_.setZero();
  P_.setIdentity();

  weights_(0) = lambda_ / (lambda_ + n_aug_);

  for (int i =1; i<2*n_aug_+1; ++i)
  {
    weights_(i) = 0.5 / (lambda_ + n_aug_); 
  }

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if (!is_initialized_) {
    // initialize the filter with measurement values and don't predict
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) 
    {     
      x_(0) = meas_package.raw_measurements_(0); // position
      x_(1) = meas_package.raw_measurements_(1); // velocity
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;  
      
      cout <<"init LASER"<<x_<< endl;
    } 
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) 
    {
      // convert coordinates from polar to cartesian
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double drho = meas_package.raw_measurements_(2);
      
      double dx = drho * cos(phi);
      double dy = drho * sin(phi);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = sqrt(pow(dx,2) + pow(dy,2));  
      x_(3) = 0;
      x_(4) = 0;   

      cout <<"init RADAR"<<x_<< endl;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

  } 
  else 
  {
    // predict and update
    // get the current and elapsed times
    double dt = (meas_package.timestamp_ - time_us_) / 1e6;
    time_us_ = meas_package.timestamp_;

    // predict using motion model
    Prediction(dt);

    // updated based on measurement type
    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) 
    { 
      UpdateLidar(meas_package);
    } 
    else if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) 
    {
      UpdateRadar(meas_package);
    }
    cout << "timestamp : " << time_us_ << endl;
    cout << "State X : " << x_.transpose() << endl;
  }
  //std::cout<<"Process Measurement End"<<std::endl;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // create augmented matrices
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  // populate augmented matrices
  x_aug.setZero();
  x_aug.head(n_x_) = x_;
  cout<<"x aug : " << x_.transpose() <<endl;

  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  //cout<<"P aug"<<P_aug<<endl;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  cout<<"Xsig_aug_ : "<< Xsig_aug <<endl;

  // transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_ + 1; ++i) 
  {
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double dyaw     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
      double px_p, py_p;
    if (fabs(dyaw) >  0.001)
    {
      px_p = p_x + v/dyaw * (sin(yaw + dyaw*delta_t) - sin(yaw)) + 0.5*nu_a*pow(delta_t,2) * cos(yaw);;
      py_p = p_y + v/dyaw * (cos(yaw) - cos(yaw + dyaw*delta_t)) + 0.5*nu_a*pow(delta_t,2) * sin(yaw);;
    } 
    else 
    {
      px_p = p_x + v*delta_t*cos(yaw) + 0.5*nu_a*pow(delta_t,2) * cos(yaw);;
      py_p = p_y + v*delta_t*sin(yaw) + 0.5*nu_a*pow(delta_t,2) * sin(yaw);;
      //std::cout <<"yawd < 0.001" << std::endl;
    }

    double v_p = v + nu_a*delta_t;
    double yaw_p = yaw + dyaw*delta_t + 0.5*nu_yawdd*pow(delta_t,2);;
    double yawd_p = dyaw + nu_yawdd*delta_t;;
    
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

  }
  cout<<"Xsig_pred : " << Xsig_pred_ <<endl;

  // predicted state mean
  x_.setZero();
  x_ = Xsig_pred_ * weights_;
  
  // predicted state covariance
  //std::cout<<"predicted state covariance"<<std::endl;

  P_.setZero();
  for (int i = 0; i < 2*n_aug_+1; ++i) 
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) >  M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
  //cout<<"P_"<<P_<<std::endl;
  //cout<<"x_ : " << x_.transpose() <<endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  
  // radar incoming measurement
  VectorXd z = VectorXd(n_z_laser_);
  z = meas_package.raw_measurements_;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);
  Zsig.setZero();
  
  // mean predicted measurement
  VectorXd z_p = VectorXd(n_z_laser_);
  z_p.setZero();
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_laser_,n_z_laser_);
  S.setZero();

  MatrixXd R = MatrixXd(n_z_laser_,n_z_laser_);
  R.setZero();
  R(0, 0) = pow(std_laspx_,2);
  R(1, 1) = pow(std_laspy_,2);
  
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);
  Tc.setZero();

  // Vector of StateDiff and MeasurmentDiff
  VectorXd x_diff  = VectorXd(n_x_);
  VectorXd z_diff = VectorXd(n_z_radar_);

  // Kalman Gain
  MatrixXd K = MatrixXd(n_x_, n_z_laser_);
  
  // transform sigma points into measurement space
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }
   // calculate mean predicted measurement
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      z_p += weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S and cross covariance matrix Tc
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
    // residual
    z_diff = Zsig.col(i) - z_p;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();

    // state difference
    x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  S += R;

  // calculate Kalman gain K;
  K = Tc*S.inverse();

  // residual
  z_diff = z - z_p;
  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();

  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // set measurement dimension, radar can measure r, phi, and r_dot
  // radar incoming measurement
  VectorXd z = VectorXd(n_z_radar_);
  z = meas_package.raw_measurements_;
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  Zsig.setZero();
  // mean predicted measurement
  VectorXd z_p = VectorXd(n_z_radar_);
  z_p.setZero();

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
  S.setZero();
  MatrixXd R = MatrixXd(n_z_radar_,n_z_radar_);
  R.setZero();
  R(0, 0) = pow(std_radr_,2);
  R(1, 1) = pow(std_radphi_,2);
  R(2, 2) = pow(std_radrd_,2);
  
  // cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  Tc.setZero();
  
  // Vector of StateDiff and MeasurmentDiff
  VectorXd x_diff  = VectorXd(n_x_);
  VectorXd z_diff = VectorXd(n_z_radar_);

  // Kalman Gain
  MatrixXd K = MatrixXd(n_x_, n_z_laser_);  

  // transform sigma points into measurement space
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    
    Zsig(0, i) = sqrt(pow(px, 2) + pow(py, 2));
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px*cos(yaw)*v + py*sin(yaw)*v) / (sqrt(pow(px, 2) + pow(py, 2)));
  }

  // calculate mean predicted measurement
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
      z_p += weights_(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S and cross covariance matrix Tc
  for (int i=0; i<2 * n_aug_ + 1; i++)
  {
    // residual
    z_diff = Zsig.col(i) - z_p;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  S += R;
  // calculate Kalman gain K;
  K = Tc*S.inverse();

  // residual
  z_diff = z - z_p;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K*(z_diff);
  P_ = P_ - K*S*K.transpose();

  while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
  while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
}
