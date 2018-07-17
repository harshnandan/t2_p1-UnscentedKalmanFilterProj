#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  previous_timestamp_ = 0;

  n_x_ = 5;
  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  is_initialized_ = false;
  
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_, 0,
	  		  0, std_laspy_*std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
			  0, std_radphi_*std_radphi_,0,
			  0, 0, std_radrd_*std_radrd_;

  NIS_ = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

   if(is_initialized_==false){
	  x_ << 0, 0, 0, 0, 0;
	  P_ <<    0.1500,    0.0000,    0.0000,    0.0000,    0.0000,
	           0.0000,    0.1500,    0.0000,    0.0000,    0.0000,
			   0.0000,    0.0000,    1.0000,    0.0000,    0.0000,
	           0.0000,    0.0000,    0.0000,    1.0000,    0.0000,
	           0.0000,    0.0000,    0.0000,    0.0000,    1.0000;

	   if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
		  float rho = meas_package.raw_measurements_(0);
		  float theta = meas_package.raw_measurements_(1);
		  x_ << rho*cos(theta), rho*sin(theta), 0, 0, 0;
	   }else{
		   float px = meas_package.raw_measurements_(0);
		   float py = meas_package.raw_measurements_(1);

		   x_ << px, py, 0, 0, 0;
	   }
	   is_initialized_ = true;
	   previous_timestamp_ = meas_package.timestamp_  ;
	   return;
   }

   double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000;
   previous_timestamp_ = meas_package.timestamp_;

   Prediction(dt);
   if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
	   if(use_radar_==true){
		   UpdateRadar(meas_package);
	   }
   }else{
	   if(use_laser_==true){
		   UpdateLidar(meas_package);
	   }
   }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++){
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // Predict mean state and process covariance matrix
  x_.fill(0.0);
  P_.fill(0.0);
  float wi;

    for(int i=0; i<2*n_aug_+1; i++){
     if (i==0){
         wi = lambda_/(lambda_ + n_aug_);
     }else{
         wi = 1/(2*(lambda_ + n_aug_));
      }
      x_ = x_ + wi * Xsig_pred_.col(i);
    }

 for(int i=0; i<2*n_aug_+1; i++){
     if (i==0){
         wi = lambda_/(lambda_ + n_aug_);
     }else{
         wi = 1/(2*(lambda_ + n_aug_));
      }
      MatrixXd diff = (Xsig_pred_.col(i) - x_);
      MatrixXd diffT = diff.transpose();
      P_ = P_ + wi * diff * diffT ;
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	  VectorXd z_pred = VectorXd(2);
	  VectorXd z = VectorXd(2);
	  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
	  z_pred.fill(0.0);

	  MatrixXd H = MatrixXd(2,5);
	  H << 1, 0, 0, 0, 0,
		   0, 1, 0, 0, 0;
	  z_pred = H * x_;

	  VectorXd y = z - z_pred;

	  MatrixXd Ht = H.transpose();
	  MatrixXd S = H * P_ * Ht + R_lidar_;

	  MatrixXd PHt = P_ * Ht;

	  MatrixXd Si = S.inverse();
	  MatrixXd K = PHt * Si;

	//new estimate
	  x_ = x_ + (K * y);
	  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
	  P_ = (I - K * H) * P_;

	  // NIS
	  NIS_ = y.transpose() * S.inverse() * y;
	  cout << NIS_ << "\n";
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	int n_z = 3;
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	VectorXd z_pred = VectorXd(3);
	z_pred.fill(0.0);

	VectorXd wi = VectorXd(2 * n_aug_ + 1);

	for(int i=0; i<2*n_aug_+1; i++){
		if (i==0){
			wi(i) = lambda_/(lambda_ + n_aug_);
		}else{
			wi(i) = 0.5/(lambda_ + n_aug_);
		}

		float px   = Xsig_pred_(0, i);
		float py   = Xsig_pred_(1, i);
		float v    = Xsig_pred_(2, i);
		float yaw  = Xsig_pred_(3, i);
		//float yawd = Xsig_pred_(4, i);

		Zsig(0, i) = sqrt(px*px + py*py);
		Zsig(1, i) = atan2(py, px);
		Zsig(2, i) = (px*cos(yaw)*v+py*sin(yaw)*v)/sqrt(px*px + py*py);

		z_pred = z_pred + wi(i) * Zsig.col(i);
	}

	MatrixXd S = MatrixXd(3, 3);
	S.fill(0.0);
	for(int i=0; i<2*n_aug_+1; i++){
		MatrixXd diff = Zsig.col(i) - z_pred;
		MatrixXd diffT = diff.transpose();
		S = S + wi(i) * (diff * diffT);
	}

	S = S + R_radar_;

	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

	    //residual
	    VectorXd z_diff = Zsig.col(i) - z_pred;
		  //angle normalization
		  while (z_diff(1)> M_PI) {
			  z_diff(1)-=2.*M_PI;
		  }
		  while (z_diff(1)<-M_PI) {
			  z_diff(1)+=2.*M_PI;
		  }

	    // state difference
	    VectorXd x_diff = Xsig_pred_.col(i) - x_;
	    //angle normalization
		  while (x_diff(3)> M_PI) {
			  x_diff(3)-=2.*M_PI;
		  }
		  while (x_diff(3)<-M_PI) {
			  x_diff(3)+=2.*M_PI;
		  }
	    //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
	    //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

	    Tc = Tc + wi(i) * x_diff * z_diff.transpose();
	  }

	  //Kalman gain K;
	  MatrixXd K = Tc * S.inverse();

	  //residual
	  VectorXd z = VectorXd(3);
	  z << meas_package.raw_measurements_(0),
			  meas_package.raw_measurements_(1),
			  meas_package.raw_measurements_(2);
	  VectorXd y = z - z_pred;

	  //angle normalization
	  while (y(1)> M_PI) {
		  y(1)-=2.*M_PI;
		  //cout << "z_diff(1): " << z_diff(1) << "\n";
	  }
	  while (y(1)<-M_PI) {
		  y(1)+=2.*M_PI;
		  //cout << "z_diff(1): " << z_diff(1) << "\n";
	  }

	  //update state mean and covariance matrix
	  x_ = x_ + K * y;
	  P_ = P_ - K*S*K.transpose();

	  // NIS
	  NIS_ = y.transpose() * S.inverse() * y;
	  cout << NIS_ << "\n";
}

