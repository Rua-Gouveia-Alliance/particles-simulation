#pragma once

#include "Particle.hpp"
#include <vector>

class CellInfo {
public:
  double x, y;
  bool counted = false;
  long id, weight = 0, rank = -1;
  std::vector<Particle> particles;

  CellInfo(long id, double x, double y) : id(id), x(x), y(y) {}

  void add_particle(Particle &p) {
    particles.push_back(p);
    ++weight;
  }
};
