#include "UKF.h"

#include <iostream>
#include <stdexcept>

UKF::UKF():
  m_x(5),
  m_firstUpdate(true),
  m_lastTimestamp(0), m_P(5, 5), m_Q(2, 2),
  m_RRadar(3, 3), m_RLidar(2, 2),
  m_nis(0) {

  m_x.fill(0);

  // Initialize P
  // Unknown position and speed
  for (int i = 0; i < 5; i++)
    m_P(i, i) = 1;

  // Initialize process noise sigma
  m_stdA = 6;
  m_stdYawdd = 6;

  m_Q << (m_stdA * m_stdA) , 0,
      0, (m_stdYawdd * m_stdYawdd);

  // Initialize R matrices
  m_RRadar(0, 0) = m_RRadar(2, 2) = 0.09;
  m_RRadar(1, 1) = 0.0009;
  m_RLidar(0, 0) = m_RLidar(1, 1) = 0.0225;

  // Initialize temp variables
  m_XSigmaAug = Eigen::MatrixXd(m_nAug, 2 * m_nAug + 1);
  m_XSigmaPred = Eigen::MatrixXd(m_nx, 2 * m_nAug + 1);
  m_xPred = Eigen::VectorXd(m_nx);
  m_PPred = Eigen::MatrixXd(m_nx, m_nx);
  m_weights = Eigen::VectorXd(2 * m_nAug + 1);
  m_weights[0] = m_lambda / (m_lambda + m_nAug);
  for (int i = 1; i < 2 * m_nAug + 1; i++)
    m_weights[i] = 0.5 / (m_lambda + m_nAug);
}

void UKF::generateSigmaPoints() {
  Eigen::VectorXd xAug(m_nAug);
  xAug.fill(0);
  xAug.head(5) = m_x;
  Eigen::MatrixXd PAug(m_nAug, m_nAug);
  PAug.topLeftCorner(5, 5) = m_P;
  PAug.bottomRightCorner(2, 2) = m_Q;
  m_XSigmaAug.col(0) = xAug;

  Eigen::MatrixXd A = PAug.llt().matrixL();

  for (int i = 0; i < m_nAug; i++) {
    float mul = sqrt(m_lambda + m_nAug);
    m_XSigmaAug.col(i + 1) = xAug + mul * PAug.col(i);
    m_XSigmaAug.col(i + 1 + m_nAug) = xAug - mul * PAug.col(i);
  }
}

void UKF::predictSigmaPoints(double delta) {
  double deltaSq = delta * delta;
  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    float px = m_XSigmaAug(0, i);
    float py = m_XSigmaAug(1, i);
    float v = m_XSigmaAug(2, i);
    float psi = m_XSigmaAug(3, i);
    float psiDot = m_XSigmaAug(4, i);
    float nuA = m_XSigmaAug(5, i);
    float nuPsiDotDot = m_XSigmaAug(6, i);

    Eigen::VectorXd xPred(m_nx);
    xPred << px, py, v, psi, psiDot;
    Eigen::VectorXd processDelta(m_nx);

    if (fabs(psiDot) >= 1e-5) {
      processDelta << (v/psiDot) * (sin(psi + psiDot * delta) - sin(psi)),
                      (v/psiDot) * (cos(psi) - cos(psi + psiDot * delta)),
                      0, psiDot * delta, 0;
    } else {
      processDelta << v * cos(psi) * delta, v * sin(psi) * delta,
                      0, 0, 0;
    }
 
    Eigen::VectorXd noiseDelta(m_nx);
    noiseDelta << 1/2.0 * deltaSq * cos(psi) * nuA,
                  1/2.0 * deltaSq * sin(psi) * nuA,
                  delta * nuA,
                  1/2.0 * deltaSq * nuPsiDotDot,
                  delta * nuPsiDotDot;
 
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
  
  //measurement covariance matrix S
  m_S = Eigen::MatrixXd(3, 3);

  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd ev = m_XSigmaPred.col(i);
    double rho = sqrt(ev[0] * ev[0] + ev[1] * ev[1]);
    double phi = fabs(ev[0]) >= 1e-5 ? atan2(ev[1], ev[0]) : 0;
    double rhoDot = fabs(rho) < 1e-5 ? 0 : ev[2] / rho * (ev[0] * cos(ev[3]) + ev[1] * sin(ev[3]));
    normalizeAngle(phi);
    m_ZSigma(0, i) = rho;
    m_ZSigma(1, i) = phi;
    m_ZSigma(2, i) = rhoDot;
    m_zPred += m_weights[i] * m_ZSigma.col(i);
  }
   
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
  
  //measurement covariance matrix S
  m_S = Eigen::MatrixXd(2, 2);

  for (int i = 0; i < 2 * m_nAug + 1; i++) {
    Eigen::VectorXd ev = m_XSigmaPred.col(i);
    m_ZSigma(0, i) = m_XSigmaPred(0, i);
    m_ZSigma(1, i) = m_XSigmaPred(1, i);
    m_zPred += m_weights[i] * m_ZSigma.col(i);
  }
   
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

  float dt = (update.getTimestamp() - m_lastTimestamp) / 1000000.0;
  m_lastTimestamp = update.getTimestamp();
  
  predictState(dt);
  predictLidarMeasurement();
  calculateGain();

  // Measure Update
  Eigen::VectorXd z(2);
  z(0) = update.getX();
  z(1) = update.getY();

  updateState(z);  
}

void UKF::processRadarUpdate(const RadarUpdate& update) {
  if (m_firstUpdate) {
    m_firstUpdate = false;
    float rho = update.getRho();
    float phi = update.getPhi();
    m_x(0) = rho * cos(phi); // px
    m_x(1) = rho * sin(phi); // py
    m_lastTimestamp = update.getTimestamp();
    return;
  }

  float dt = (update.getTimestamp() - m_lastTimestamp) / 1000000.0;
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
}
