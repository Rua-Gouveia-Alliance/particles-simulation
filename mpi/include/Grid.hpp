#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "Cell.hpp"
#include "Mass.hpp"
#include "Particle.hpp"
#include <unordered_map>
#include <vector>

class Grid {
private:
  double _side;
  long _ncside;
  Particle &_first_particle;
  Particle _default_first_particle = Particle();
  std::unordered_map<long, Mass> _adjacent_masses;
  std::unordered_map<long, std::vector<Particle>> _new_particles;

  long _get_particle_cell(Particle &p);
  void _add_particle_to_cell(Particle &p);

public:
  std::unordered_map<long, Cell> _cells;

  Grid(double side, long ncside, std::vector<int> depedents);

  void add_cell(Cell &c);
  std::vector<long> get_adjacent_cells(long ci);
  void update_cells();
  void print_cells() const;
  Particle get_first_particle() const;
  long get_collisions();
};

#endif
