#include "Particle.hpp"

#include <iostream>

Particle::Particle() : x(0), y(0), vx(0), vy(0), m(0), first_particle(false) {}

Particle::Particle(double m)
    : x(0), y(0), vx(0), vy(0), m(m), first_particle(false) {}

Particle::Particle(double x, double y, double vx, double vy, double m)
    : x(x), y(y), vx(vx), vy(vy), m(m), first_particle(false) {}

void Particle::print_info() const {
  std::cout << "x: " << x << "\ty: " << y << "\tvx: " << vx << "\tvy: " << vy
            << "\tm: " << m << std::endl;
}
