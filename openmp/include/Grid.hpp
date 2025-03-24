#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "Cell.hpp"
#include "Particle.hpp"
#include <vector>

class Grid {
private:
  double _side;
  long _ncside;
  Particle &_first_particle;
  // TODO: Possibly not the best workaround
  Particle _default_first_particle = Particle();

  void _add_particle_to_cell(Particle &p);

public:
  std::vector<Cell> _cells;

  Grid(double side, long ncside);

  void add_cell(Cell &c);
  std::vector<long> get_adjacent_cells(long ci);
  void update_cells();
  void print_cells() const;
  Particle get_first_particle() const;
  long get_collisions();
};

#endif
