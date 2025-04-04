#pragma once

#include "CellInfo.hpp"
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

  PartialCell(CellInfo info, double side)
      : _id(info.id), owner(info.rank), x(info.x), y(info.y), side(side),
        mass({info.id, -1, -1, 0}){};

  void add_particle(Particle &p) { particles.push_back(p); }
  void clear() { particles.clear(); }

  int id() const { return _id; }
};
