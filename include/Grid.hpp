#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "Cell.hpp"
#include <vector>

class Grid {
private:
  double _side;
  long _ncside;

  void _add_particle_to_cell(Particle &p);

public:
  std::vector<Cell> _cells;

  Grid(long side, long ncside);

  void add_cell(Cell &c);
  std::vector<long> get_adjacent_cells(long ci);
  void update_cells();
};

#endif
