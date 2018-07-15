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
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	cout << "estimations.size() " << estimations.size() << "\n";
	cout << "ground_truth.size() " << ground_truth.size() << "\n";
	if (estimations.size() != ground_truth.size()){
		cout << "Tools:RMSE - Estimation and Ground Truth size do not match\n";
	}

	int n = ground_truth.size();

	for(int i=0; i<n; i++){
		VectorXd err = estimations[i] - ground_truth[i];
		err = err.array() * err.array();
		rmse += err;
	}
	rmse = rmse/n;
	rmse = rmse.array().sqrt();
	return rmse;
}
