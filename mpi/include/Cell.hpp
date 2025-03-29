#pragma once

#include "Mass.hpp"
#include "PartialCell.hpp"
#include "Particle.hpp"
#include <unordered_map>
#include <vector>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

class Cell {
private:
  long _id;
  std::vector<Particle> _particles;
  std::vector<Particle> _temp_particles;

  void _check_collisions();

public:
  Mass mass;
  long collisions = 0;
  double x, y, side;
  std::vector<int> adjacent_ranks;

  Cell(long id, std::vector<int> &adjacent_ranks, double x, double y,
       double side)
      : _id(id), adjacent_ranks(adjacent_ranks), x(x), y(y), side(side),
        mass({id, -1, -1, 0}){};

  void update_mass();
  void add_particle(Particle &p);
  std::vector<Particle>
  update_particles(const std::unordered_map<int, PartialCell> &adjacent_cells);
  void finish_update();
  void print_particles() const;

  long id() const { return _id; }
  const std::vector<Particle> &get_particles() const;
};
