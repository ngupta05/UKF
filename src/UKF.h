#pragma once
#include "Update.h"
#include "Eigen/Dense"

#include <stdint.h>
#include <vector>

class UKF {
  bool m_firstUpdate;
  uint64_t m_lastTimestamp;
  Eigen::VectorXd m_x;
  Eigen::MatrixXd m_P;
  Eigen::MatrixXd m_Q;
  Eigen::MatrixXd m_RRadar;
  Eigen::MatrixXd m_RLidar;
  int m_stdA;
  int m_stdYawdd;
  float m_nisLidar;
  float m_nisRadar;

  // Temp variables for calculations
  const int m_nx = 5;
  const int m_nAug = 7;
  const float m_lambda = 3 - m_nAug;
  const int m_nzRadar = 3;
  const int m_nzLidar = 2;
  Eigen::VectorXd m_weights;
  Eigen::MatrixXd m_XSigmaAug, m_XSigmaPred;
  Eigen::VectorXd m_xPred;
  Eigen::MatrixXd m_PPred;
  Eigen::MatrixXd m_ZSigma;
  Eigen::MatrixXd m_zPred;
  Eigen::MatrixXd m_S;
  Eigen::MatrixXd m_K;
  Eigen::MatrixXd m_T;

  public:
  UKF();
  void processLidarUpdate(const LidarUpdate& update);
  void processRadarUpdate(const RadarUpdate& update);
  Eigen::VectorXd getEstimate() const { return m_x; }
  float getNISLidar() const { return m_nisLidar; }
  float getNISRadar() const { return m_nisRadar; }

  private:
  template<typename T>
  void normalizeAngle(T& phi) {
    while (phi > M_PI) phi -= 2 * M_PI;
    while (phi < -M_PI) phi += 2 * M_PI;
  }

  void generateSigmaPoints();
  void predictSigmaPoints(double delta);
  void predictMeanAndCov();
  void predictState(double delta);
  void predictRadarMeasurement();
  void predictLidarMeasurement();
  void calculateGain();
  void updateState(const Eigen::VectorXd& z);
};
