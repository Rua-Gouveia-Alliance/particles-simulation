#pragma once
#include <iostream>

class Particle {
private:
public:
  double x, y, vx, vy, m;
  bool first_particle;

  Particle() : x(0), y(0), vx(0), vy(0), m(0), first_particle(false) {}
  Particle(double m) : x(0), y(0), vx(0), vy(0), m(m), first_particle(false) {}
  Particle(double x, double y, double vx, double vy, double m)
      : x(x), y(y), vx(vx), vy(vy), m(m), first_particle(false) {}

  void print_info() const {
    std::cout << "x: " << x << "\ty: " << y << "\tvx: " << vx << "\tvy: " << vy
              << "\tm: " << m << std::endl;
  }
};
