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
    
    /*
     * Jacobian marix, init it as part of
     * FusionEKF obj and set it on every
     * radar measurement
     */
    Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

  //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;
 
    /*
     * Set the process and measurement noises
     */
    // Laser measurement function
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    //State vector
    ekf_.x_ = VectorXd(4);


    //State covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
   
    /*
     * Initialize the covariance  with
     * high uncertainity for vx and vy
     *
     * tried initialized with low uncertainity for
     * vx and vy for radar data but wasn't make
     * significant difference
     */
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0, // vx
               0, 0, 0, 1000; // vy


    //Initial State transition Matrix.
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    ekf_.Q_ = MatrixXd(4, 4);
   
    // Set acceleration noise components
    noise_ax = 9;
    noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

/*
 * Calculate Process Coveriance Matrix Q
 */
void
FusionEKF::UpdateCovarianceQ (KalmanFilter &kf_, float dt)
{
    float dt_2, dt_3, dt_4;
    
    dt_2 = dt * dt;
    dt_3 = dt_2 * dt;
    dt_4 = dt_3 * dt;

    /*
     * set the Coveriance Matrix Q into KF.
     */
    kf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
               0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
               dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
               0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
}
void 
FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
    float dt; 
    float ro, phi, ro_dot; /* polar coordinate */
    float px, py, vx, vy;  /* cartesian coordinate */
    Tools tool;            /* for calculation jacobian m and rmse */ 
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    
    if (!is_initialized_) {

        /*
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        cout << "EKF: " << endl;
       
        previous_timestamp_ = measurement_pack.timestamp_;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
          
            /* first measurement */

            /*
             * Radar data isin polar coordinate
             * Convert to cartesian coordinates and initialize state.
             *
             * Map initial motion state from polar to cartesian  coordinate
             */
            ro = measurement_pack.raw_measurements_[0];
            phi = measurement_pack.raw_measurements_[1];
            ro_dot = measurement_pack.raw_measurements_[2];

            px = ro*cos(phi);
            py = ro*sin(phi);
            vx = ro_dot*cos(phi);
            vy = ro_dot*sin(phi);
            ekf_.x_ << px, py, vx, vy;

            /*
             * lets initilize when a radar data is available
             * to have an initial vx and vy 
             */
            is_initialized_ = true;

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
          
            /*
             * Lidar is in cartesian axis and can not measure velocity, so just
             * put vx, vy = 0, 0
             * Initialize motion state.
             */
            px  = measurement_pack.raw_measurements_[0];
            py = measurement_pack.raw_measurements_[1]; 
            ekf_.x_ << px, py, 0, 0;
       }

        // done initializing, no need to predict or update
        //is_initialized_ = true;
        return;
    }

    /*
     * compute the time elapsed between the current and previous measurements
     * dt - expressed in seconds
     */
    dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    /*
     * Update the state transition matrix F according to the new elapsed time.
     * - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    /*
     * update F matrix to integrate the delta time
     */
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    /*
     * Update the process noise covariance matrix.
     */
     UpdateCovarianceQ(ekf_, dt);
    
    /*
     * Kalman filter predict
     */
     ekf_.Predict();


  /*****************************************************************************
   *  Update
   ****************************************************************************/

    /*
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        
        /*
         * Measurement from Radar, so calculate the new jacobian matrix Hj,
         * use the non-linear measurement function (contrast with linear
         * measurement function H for laser) to project the predictate
         * state
         */
        cout << "RADAR" << endl;

        /*
         * so Basically using Jacobian matrix we are mapping the predicted
         * state 4 dimension (px, py', vx', vy') into Radar measurement space
         * of 3 dimension (ro, phi, ro_dot) and then update the measurement
         */
        tool.CalculateJacobian(ekf_.x_, Hj_);
        ekf_.H_ = Hj_;
        ekf_.R_ =  R_radar_;
        
        /*
         * Radar updates
         */
         ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    } else {
        
        /*
         * Laser updates
         */
        cout << "LASER" << endl;

        ekf_.R_ = R_laser_;
        ekf_.H_ =  H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
