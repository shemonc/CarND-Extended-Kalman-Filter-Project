#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  
    /*
     * predict the state. Same for both lader and Radar in EKF
     */
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  
	/*
     * update the state by using Kalman Filter equations
  	 */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	/*
     * new estimate
     */
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
    float ro, phi, ro_dot;
    float px, py, vx, vy;
    
    /*
     * update the state by using Extended Kalman Filter equations
  	 */

    /*
     * First convert predictate motion state x(t+1) at time t from
     * cartesian to polar coordinate
     */
    px = x_(0);
    py = x_(1);
    vx = x_(2);
    vy = x_(3);
    ro = sqrt(px*px + py*py);
    phi = atan2(py,px);

    /*
     * if range is close to zero then set the range rate to zero.
     * no change 
     */
    if (fabs(ro) < 0.0001) {
        ro_dot = 0;
    } else {
        ro_dot = (px*vx + py*vy) / ro;
    }

    /*
     * error measurement
     * difference between measurement directly from radar and predictate state
     * both in polar axis
     */
	VectorXd z_pred(3);
	z_pred << ro, phi, ro_dot;
    VectorXd y = z - z_pred;

    /*
     * phi adjustment
     * atan2() returns values between -pi and pi, When calculating phi in 
     * y = z - h(x) for radar measurements, the resulting angle phi in the y
     * vector should be adjusted so that it is between -pi and pi. The Kalman
     * filter is expecting small angle values between the range -pi and pi.
     * Add 2π or subtract 2π until the angle is within the desired range.
     */
    if (y(1) > M_PI) y(1) -= 2*M_PI;
    if (y(1) < -M_PI) y(1) += 2*M_PI;

    /*
     * H_ is Jacobian matrix
     */
    MatrixXd Ht = H_.transpose();
	
    /*
     * error matrix
     */
    MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
