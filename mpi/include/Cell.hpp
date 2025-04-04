#pragma once

#include "CellInfo.hpp"
#include "Mass.hpp"
#include "PartialCell.hpp"
#include "Particle.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

class Cell {
private:
  long _id;
  std::vector<Particle> _temp_particles;

public:
  Mass mass;
  double x, y, side;
  long collisions = 0;
  std::vector<Particle> particles;
  std::unordered_set<int> adjacent_ranks;

  Cell(CellInfo info, std::unordered_set<int> &adjacent_ranks, double side)
      : _id(info.id), x(info.x), y(info.y),
        adjacent_ranks(std::move(adjacent_ranks)), side(side),
        mass({info.id, -1, -1, 0}), particles(std::move(info.particles)){};

  void update_mass();
  void add_particle(Particle &p) { _temp_particles.push_back(p); }
  void add_updated_particle(Particle &p) { particles.push_back(p); }
  std::vector<Particle>
  update_particles(const std::vector<Mass> &adjacent_masses);
  void finish_update();
  void check_collisions();
  void print_particles() const;

  long id() const { return _id; }
};
