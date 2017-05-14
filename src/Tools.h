#pragma once
#include "Eigen/Dense"
#include <vector>

class Tools {
  public:
    static Eigen::VectorXd calculateRMSE(std::vector<Eigen::VectorXd> estimates,
        std::vector<Eigen::VectorXd> truth);
};
