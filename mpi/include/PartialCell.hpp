#pragma once

#include "Mass.hpp"
#include "Particle.hpp"
#include <vector>

class PartialCell {
private:
  long _id;
  int _owner;
  std::vector<Particle> _particles;

public:
  Mass mass;
  double x, y, side;

  PartialCell(long id, int owner, double x, double y, double side)
      : _id(id), _owner(owner), x(x), y(y), side(side), mass({id, -1, -1, 0}){};

  void add_particle(Particle &p) { _particles.push_back(p); }
  void clear() { _particles.clear(); }

  int id() const { return _id; }
  const std::vector<Particle> &get_particles() const { return _particles; }
};
