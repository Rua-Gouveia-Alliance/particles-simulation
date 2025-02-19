#include "Particle.hpp"

Particle::Particle() : _x(0), _y(0), _vx(0), _vy(0), _m(0) {}

Particle::Particle(double x, double y, double vx, double vy, double m)
    : _x(x), _y(y), _vx(vx), _vy(vy), _m(m) {}
