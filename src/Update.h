#pragma once

#include <stdint.h>

class RadarUpdate {
  float m_rho;
  float m_phi;
  float m_rhodot;
  uint64_t m_timestamp;

  public:
  RadarUpdate(float rho, float phi, float rhodot, uint64_t timestamp):
    m_rho(rho), m_phi(phi), m_rhodot(rhodot), m_timestamp(timestamp) {
  }

  float getRho() const {
    return m_rho;
  }

  float getPhi() const {
    return m_phi;
  }

  float getRhoDot() const {
    return m_rhodot;
  }

  uint64_t getTimestamp() const {
    return m_timestamp;
  }
};

class LidarUpdate {
  float m_x;
  float m_y;
  uint64_t m_timestamp;

  public:
  LidarUpdate(float x, float y, uint64_t timestamp):
    m_x(x), m_y(y), m_timestamp(timestamp) {
  }

  float getX() const {
    return m_x;
  }

  float getY() const {
    return m_y;
  }

  uint64_t getTimestamp() const {
    return m_timestamp;
  }
};
