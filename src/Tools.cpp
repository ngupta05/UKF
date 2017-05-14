#include "Tools.h"

#include <iostream>

Eigen::VectorXd Tools::calculateRMSE(std::vector<Eigen::VectorXd> estimates,
    std::vector<Eigen::VectorXd> truth) {

  if (0 == estimates.size() || estimates.size() != truth.size()) {
    std::cout << "Invalid input" << std::endl;
    return Eigen::VectorXd();
  }
  unsigned int stateSize = estimates[0].size();

  if (0 == stateSize)
    return Eigen::VectorXd();

  Eigen::VectorXd rmse(stateSize);
  rmse.fill(0);
  unsigned int size = estimates.size();
  for (unsigned int i = 0; i < size; i++) {
    Eigen::VectorXd diff = estimates[i] - truth[i];
    diff = diff.array() * diff.array();
    rmse += diff;
  }

  rmse /= size;

  rmse = rmse.array().sqrt();
  return rmse;
}
