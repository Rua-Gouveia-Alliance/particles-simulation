#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "Cell.hpp"
#include "Particle.hpp"
#include <vector>

class Grid {
private:
  double _side;
  long _ncside;

  void _add_particle_to_cell(Particle &p);

public:
  std::vector<Cell> _cells;

  Grid(double side, long ncside);

  void add_cell(Cell &c);
  std::vector<long> get_adjacent_cells(long ci);
  void update_cells();
  void print_cells();
  Particle get_first_particle();
  long get_collisions();
};

#endif
