#include "Particle.hpp"

#include <iostream>

Particle::Particle()
    : _x(0), _y(0), _vx(0), _vy(0), _m(0), _first_particle(false) {}

Particle::Particle(double m)
    : _x(0), _y(0), _vx(0), _vy(0), _m(m), _first_particle(false) {}

Particle::Particle(double x, double y, double vx, double vy, double m)
    : _x(x), _y(y), _vx(vx), _vy(vy), _m(m), _first_particle(false) {}

void Particle::print_info() {
  std::cout << "x: " << _x << "\ty: " << _y << "\tvx: " << _vx
            << "\tvy: " << _vy << "\tm: " << _m << std::endl;
}
