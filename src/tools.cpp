#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

/*
 * CalculateRMSE
 * 
 * Calculate root mean square error
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	
	rmse << 0,0,0,0;

	/*
	 * check the validity of the following inputs:
	 * the estimation vector size should not be zero
	 * the estimation vector size should equal ground truth vector size
	 */
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	/*
	 * accumulate squared residuals
     */
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

/*
 * CalculateJacobian
 *
 * Calculate the Jacobian matrix
 * @in, x_state predictate state 
 * @in, Hj Ref to Jacobian matrix obj to be calculated
 */
void
Tools::CalculateJacobian(const VectorXd& x_state, MatrixXd &Hj) {

    float px, py, vx, vy;
    float c1, c2, c3;

    
    /*
     * Reset/init to all zeros
     */
    Hj << 0,0,0,0,
          0,0,0,0,
          0,0,0,0;

    /* retrieve the state parameters */
    px = x_state(0);
    py = x_state(1);
    vx = x_state(2);
    vy = x_state(3);

    /* pre-compute a set of terms to avoid repeated calculation */
    c1 = px*px + py*py;
    c2 = sqrt(c1);
    c3 = (c1*c2);

    if (fabs(c1) < 0.0001) {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return;
    }

    /* compute the Jacobian matrix */
    Hj << (px/c2), (py/c2), 0, 0,
         -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
}
