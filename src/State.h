#pragma once
#include <math.h>

struct State {
  float x;
  float y;
  float vx;
  float vy;

  State() : x(0), y(0), vx(0), vy(0) {
  }

  State operator-(const State& rhs) const {
    State diff;
    diff.x = x - rhs.x;
    diff.y = y - rhs.y;
    diff.vx = vx - rhs.vx;
    diff.vy = vy - rhs.vy;
    return diff;
  }

  void square() {
    x = x * x;
    y = y * y;
    vx = vx * vx;
    vy = vy * vy;
  }

  void sqrt() {
    x = ::sqrt(x);
    y = ::sqrt(y);
    vx = ::sqrt(vx);
    vy = ::sqrt(vy);
  }

  State& operator+=(const State& rhs) {
    x += rhs.x;
    y += rhs.y;
    vx += rhs.vx;
    vy += rhs.vy;
    return *this;
  }
};
