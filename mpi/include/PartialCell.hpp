#pragma once

#include "Mass.hpp"
#include "Particle.hpp"
#include <vector>

class PartialCell {
private:
  long _id;

public:
  int owner;
  Mass mass;
  double x, y, side;
  std::vector<Particle> particles;

  PartialCell(long id, int owner, double x, double y, double side)
      : _id(id), owner(owner), x(x), y(y), side(side), mass({id, -1, -1, 0}){};

  void add_particle(Particle &p) { particles.push_back(p); }
  void clear() { particles.clear(); }

  int id() const { return _id; }
};
