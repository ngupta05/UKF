#include "UKF.h"

#include <iostream>
#include <stdexcept>

UKF::UKF():
  m_x(5),
  m_firstUpdate(true),
  m_lastTimestamp(0), m_P(5, 5), m_Q(2, 2),
  m_RRadar(3, 3), m_RLidar(2, 2),
  m_nis(0),
  m_xAug(7), m_PAug(7, 7),
  m_XSigmaAug(m_nAug, 2 * m_nAug + 1), 
  m_XSigmaPred(m_nx, 2 * m_nAug + 1),
  m_xPred(m_nx),
  m_PPred(m_nx, m_nx),
  m_weights(2 * m_nAug + 1) {
  m_x.fill(0);

  // Initialize P
  // Unknown position and speed
  m_P.fill(0);
  for (int i = 0; i < 5; i++)
    m_P(i, i) = 1;

  // Initialize process noise sigma
  m_stdA = 0.5;
  m_stdYawdd = 0.01;

  m_Q.fill(0);
  m_Q << (m_stdA * m_stdA) , 0,
      0, (m_stdYawdd * m_stdYawdd);

  // Initialize R matrices
  m_RRadar.fill(0);
  m_RRadar(0, 0) = m_RRadar(2, 2) = 0.09;
  m_RRadar(1, 1) = 0.0009;
  m_RLidar.fill(0);
  m_RLidar(0, 0) = m_RLidar(1, 1) = 0.0225;

  // Initialize temp variables
  m_weights[0] = m_lambda / (m_lambda + m_nAug);
  for (int i = 1; i < 2 * m_nAug + 1; i++)
    m_weights[i] = 0.5 / (m_lambda + m_nAug);
}

void UKF::generateSigmaPoints() {
  m_xAug.fill(0);
  m_xAug.head(5) = m_x;
  m_PAug.fill(0);
  m_PAug.topLeftCorner(5, 5) = m_P;
  m_PAug.bottomRightCorner(2, 2) = m_Q;
  m_XSigmaAug.col(0) = m_xAug;

  Eigen::MatrixXd A = m_PAug.llt().matrixL();

  double mul = sqrt(m_lambda + m_nAug);
  for (int i = 0; i < m_nAug; i++) {
    m_XSigmaAug.col(i + 1) = m_xAug + mul * A.col(i);
    m_XSigmaAug.col(i + 1 + m_nAug) = m_xAug - mul * A.col(i);
  }
}

void UKF::predictSigmaPoints(double delta) {
  double deltaSq = delta * delta;
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    double px = m_XSigmaAug(0, i);
    double py = m_XSigmaAug(1, i);
    double v = m_XSigmaAug(2, i);
    double yaw = m_XSigmaAug(3, i);
    double yawDot = m_XSigmaAug(4, i);
    double nuA = m_XSigmaAug(5, i);
    double nuYawDotDot = m_XSigmaAug(6, i);

    Eigen::VectorXd xPred(m_nx);
    xPred << px, py, v, yaw, yawDot;
    Eigen::VectorXd processDelta(m_nx);

    if (fabs(yawDot) >= 1e-5) {
      processDelta << (v/yawDot) * (sin(yaw + yawDot * delta) - sin(yaw)),
                      (v/yawDot) * (cos(yaw) - cos(yaw + yawDot * delta)),
                      0, yawDot * delta, 0;
    } else {
      processDelta << v * cos(yaw) * delta, v * sin(yaw) * delta,
                      0, 0, 0;
    }
 
    Eigen::VectorXd noiseDelta(m_nx);
    noiseDelta << 1/2.0 * deltaSq * cos(yaw) * nuA,
                  1/2.0 * deltaSq * sin(yaw) * nuA,
                  delta * nuA,
                  1/2.0 * deltaSq * nuYawDotDot,
                  delta * nuYawDotDot;
 
    xPred += (processDelta + noiseDelta);
      
    m_XSigmaPred.col(i) = xPred;
  }
}

void UKF::predictMeanAndCov() {
  m_xPred.fill(0);
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
      m_xPred += m_weights[i] * m_XSigmaPred.col(i);
  }

  m_PPred.fill(0);
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd v = m_XSigmaPred.col(i) - m_xPred;
    normalizeAngle(v(3));
    m_PPred += m_weights[i] * v * v.transpose();
  }
}

void UKF::predictState(double delta) {
  generateSigmaPoints();
  predictSigmaPoints(delta);
  predictMeanAndCov();
}

void UKF::predictRadarMeasurement() {
  m_ZSigma = Eigen::MatrixXd(3, 2 * m_nAug + 1);

  //mean predicted measurement
  m_zPred = Eigen::VectorXd(3);
  m_zPred.fill(0);
  
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd ev = m_XSigmaPred.col(i);
    double rho = sqrt(ev[0] * ev[0] + ev[1] * ev[1]);
    double phi = atan2(ev[1], ev[0]); // TODO handle 0 denom
    double rhoDot = ev[2] * (ev[0] * cos(ev[3]) + ev[1] * sin(ev[3])) / std::max(1e-5, rho);
    normalizeAngle(phi);
    m_ZSigma(0, i) = rho;
    m_ZSigma(1, i) = phi;
    m_ZSigma(2, i) = rhoDot;
    m_zPred += m_weights[i] * m_ZSigma.col(i);
  }
   
  //measurement covariance matrix S
  m_S = Eigen::MatrixXd(3, 3);
  m_S.fill(0);

  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd v = m_ZSigma.col(i) - m_zPred;
    normalizeAngle(v(1));
    m_S += m_weights[i] * v * v.transpose();
  }

  m_S += m_RRadar;
  m_T = Eigen::MatrixXd(m_nx, 3);
}

void UKF::predictLidarMeasurement() {
  m_ZSigma = Eigen::MatrixXd(2, 2 * m_nAug + 1);

  //mean predicted measurement
  m_zPred = Eigen::VectorXd(2);
  m_zPred.fill(0);
  
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd ev = m_XSigmaPred.col(i);
    m_ZSigma(0, i) = m_XSigmaPred(0, i);
    m_ZSigma(1, i) = m_XSigmaPred(1, i); // TODO this can be a matrix op
    m_zPred += m_weights[i] * m_ZSigma.col(i);
  }
   
  //measurement covariance matrix S
  m_S = Eigen::MatrixXd(2, 2);
  m_S.fill(0);

  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd v = m_ZSigma.col(i) - m_zPred;
    m_S += m_weights[i] * v * v.transpose();
  }

  m_S += m_RLidar;
  m_T = Eigen::MatrixXd(m_nx, 2);
}

void UKF::calculateGain() {
  m_T.fill(0);
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    m_T += m_weights(i) * (m_XSigmaPred.col(i) - m_xPred) * (m_ZSigma.col(i) - m_zPred).transpose();
  }

  //calculate Kalman gain K;
  m_K = m_T * m_S.inverse();
}

void UKF::updateState(const Eigen::VectorXd& z) {
  m_x = m_xPred + m_K * (z - m_zPred);
  m_P = m_PPred - m_K * m_S * m_K.transpose();
}

void UKF::calculateNIS(const Eigen::VectorXd& z) {
  m_nis = (z - m_zPred).transpose() * m_S.inverse() * (z - m_zPred);
}

void UKF::processLidarUpdate(const LidarUpdate& update) {
  if (m_firstUpdate) {
    m_firstUpdate = false;
    m_x(0) = update.getX();
    m_x(1) = update.getY();
    m_lastTimestamp = update.getTimestamp();
    return;
  }

  double dt = (update.getTimestamp() - m_lastTimestamp) / 1000000.0;
  m_lastTimestamp = update.getTimestamp();
  
  predictState(dt);
  predictLidarMeasurement();
  calculateGain();

  // Measure Update
  Eigen::VectorXd z(2);
  z(0) = update.getX();
  z(1) = update.getY();

  updateState(z);  
  calculateNIS(z);
}

void UKF::processRadarUpdate(const RadarUpdate& update) {
  if (m_firstUpdate) {
    m_firstUpdate = false;
    double rho = update.getRho();
    double phi = update.getPhi();
    m_x(0) = rho * cos(phi); // px
    m_x(1) = rho * sin(phi); // py
    m_lastTimestamp = update.getTimestamp();
    return;
  }

  double dt = (update.getTimestamp() - m_lastTimestamp) / 1000000.0;
  m_lastTimestamp = update.getTimestamp();
  predictState(dt);
  predictRadarMeasurement();
  calculateGain();

  // apply measurement
  Eigen::VectorXd z(3);
  z(0) = update.getRho();
  z(1) = update.getPhi();
  z(2) = update.getRhoDot();

  updateState(z);
  calculateNIS(z);
}
