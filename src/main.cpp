#include "UKF.h"
#include "Tools.h"
#include "Update.h"
#include "Eigen/Dense"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

int main(int argc, char** argv) {
  std::string usage = argv[0];
  usage += " input.txt";
  usage += " output.txt";
  if (3 != argc) {
    std::cout << usage << std::endl;
    return 0;
  }

  const char* input = argv[1];
  const char* output = argv[2];

  std::ifstream inputFile(input, std::ifstream::in);
  std::ofstream outputFile(output, std::ofstream::out);

  std::string line;
  std::string type;
  uint64_t timestamp;
  std::vector<Eigen::VectorXd> groundTruth;
  std::vector<Eigen::VectorXd> estimates;
  UKF ukf;

  // print column names for output file
  outputFile << "time_stamp" << "\t"
    << "px_state" << "\t"
    << "py_state" << "\t"
    << "v_state" << "\t"
    << "yaw_angle_state" << "\t"
    << "yaw_rate_state" << "\t"
    << "sensor_type" << "\t"
    << "NIS" << "\t"
    << "px_measured" << "\t"
    << "py_measured" << "\t"
    << "px_ground_truth" << "\t"
    << "py_ground_truth" << "\t"
    << "vx_ground_truth" << "\t"
    << "vy_ground_truth" << "\n";

  while (getline(inputFile, line)) {
    std::istringstream iss(line);

    iss >> type;

    float px, py;
    if (type.compare("L") == 0) {
      iss >> px;
      iss >> py;
      iss >> timestamp;
      LidarUpdate update(px, py, timestamp);
      ukf.processLidarUpdate(update);
    } else {
      float rho, phi, rhodot;
      iss >> rho;
      iss >> phi;
      iss >> rhodot;
      iss >> timestamp;
      RadarUpdate update(rho, phi, rhodot, timestamp);
      ukf.processRadarUpdate(update);
      px = rho * cos(phi);
      py = rho * sin(phi);
    }

    Eigen::VectorXd estimate = ukf.getEstimate();

    float pxGt, pyGt, vxGt, vyGt;
    iss >> pxGt >> pyGt >> vxGt >> vyGt;
    Eigen::VectorXd vec(4);
    vec << pxGt, pyGt, vxGt, vyGt;
    groundTruth.push_back(vec);

    // print timestamp
    outputFile << timestamp << "\t";
    // print px, py, v, yaw, yaw_rate
    for (int i = 0; i < 5; i++)
      outputFile << estimate(i) << "\t";
    // print sensor type and NIS
    if (type.compare("L") == 0) {
      outputFile << "lidar" << "\t";
      outputFile << ukf.getNISLidar() << "\t";
    } else {
      outputFile << "radar" << "\t";
      outputFile << ukf.getNISRadar() << "\t";
    }
    // print measured px, py
    outputFile << px << "\t" << py << "\t";
    // print ground truth
    outputFile << pxGt << "\t"
     << pyGt << "\t" 
     << vxGt << "\t" 
     << vyGt << "\n"; 

    // get cartesian estimate
    Eigen::VectorXd cartEstimate(4);
    cartEstimate << estimate(0), estimate(1),
                 estimate(2) * cos(estimate(3)),
                 estimate(2) * sin(estimate(3));
    estimates.push_back(cartEstimate);
  }

  Eigen::VectorXd rmse = Tools::calculateRMSE(estimates, groundTruth); 

  std::cout << "RMSE" << std::endl
    << rmse(0) << std::endl << rmse(1) << std::endl
    << rmse(2) << std::endl << rmse(3) << std::endl;

  // close files
  if (outputFile.is_open()) {
    outputFile.close();
  }

  if (inputFile.is_open()) {
    inputFile.close();
  }

  return 0;
}
