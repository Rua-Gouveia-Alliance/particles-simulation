#pragma once

#include "Mass.hpp"
#include "Particle.hpp"
#include <vector>

class PartialCell {
private:
  int _id; // n devia ser long?
  int _owner;
  std::vector<Particle> _particles;

public:
  Mass mass;
  double x, y, side;

  PartialCell(int id, int owner, double x, double y, double side)
      : _id(id), _owner(owner), x(x), y(y), side(side), mass(){}; //TODO retirei para compilar

  void add_particle(Particle &p) { _particles.push_back(p); }
  void clear() { _particles.clear(); }

  int id() const { return _id; }
  const std::vector<Particle> &get_particles() const { return _particles; }
};
